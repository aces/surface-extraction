/*-------------------------------------------------------------------------- 

File: skeletonization.c

This program skeletonizes a minc volume, based on all non minimum value voxels 
in the input volume.  The program is an implementation of the 3D 2-subfield medial
surface thinning algorithm presented in the following paper:

A medial-surface oriented 3-d two-subfield thinning algorithm
Cerng-Min Ma, Shu-Yen Wan
Pattern Recognition Letters 22 (2001) 1439-1446

Created: April 11, 2003 (Vivek Singh)

Modified: January 21, 2004 (Vivek Singh)
Replaced all MACRO volume access calls with function calls.

-------------------------------------------------------------------------- */

#include <volume_io.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time_stamp.h>
#include "skeletonization.h"



int main(int argc, char* argv[]) {
  
  char* history;
  Volume vol1, vol2,vtemp;
  Real voxel_value;
  int i,j,k;
  int sizes[3];
  char** dnames;
  
  if (( argc != 3) || (!strcmp(argv[1],"-usage")) || (!strcmp(argv[1],"-help")))  {
    
    if ((argc > 1) &&(!strcmp(argv[1],"-help")))
      printf("Skeletonizes a binary thresholded minc volume\n");

    printf("usage: %s <input.mnc> <output.mnc>\n\n",argv[0]); 
    return ( 1 );
  }
  
  //input volume to be skeletonized
  if( input_volume( argv[1], 3, NULL, NC_UNSPECIFIED, FALSE,
		    0.0,0.0, TRUE, &vol1,
		    (minc_input_options *) NULL ) != OK )
    return( 1 );
  
  history = time_stamp(argc, argv);
  
  get_volume_sizes( vol1, sizes );
  
  vol2 = copy_volume_definition(vol1, NC_UNSPECIFIED, TRUE,0.0, 0.0);
  
  set_volume_real_range(vol2,0,1);

  //initialize voxels to zero in the output volume
  for( k = 0;  k < sizes[0];  k++ ) {  
    for( j = 0;  j < sizes[1];  j++ ) {
      for( i = 0;  i < sizes[2];  i++ ) {
	set_volume_real_value(vol2,k,j,i,0,0,VALUE_0);
      }
    }
  }

  //skeletonize vol1 to vol2
  reducer(vol1, vol2);  

  //output skeleton
  output_modified_volume( argv[2], NC_UNSPECIFIED,
	 		FALSE, 0.0, 0.0, vol2, argv[1],
			history,
			(minc_output_options *) NULL );  
  printf("\n");
  return 0;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : reducer
@INPUT      : vol1
@OUTPUT     : vol2
@RETURNS    : 
@DESCRIPTION: Reduces vol1 to a skeleton outputted in vol2
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */

void reducer(Volume vol1, Volume vol2) {

  int change = 1;
  int i;
  int n = 0;

  //keep reducing the input volume until no voxels are removed
  while(change) {

    int cont = 0;
    int tempchange[MAX_SUB] = { 0, 0 };

    printf("ITERATION %d\n",++n);

    //the algorithm uses a two subfield approach, therefore, 
    //the volumes are reduced alternatively
    tempchange[0] = sub_reducer(vol1,vol2,0);
    tempchange[1] = sub_reducer(vol2,vol1,1);
    
    for(i = 0; i < MAX_SUB; i++) {
      if (tempchange[i] == 1) {
	cont = 1;
	break;
      }
    }
    if (!cont)
      break;

  }
    
  //final step where the ends of the skeleton are removed
  reduce_ends(vol1,vol2);

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : sub_reducer
@INPUT      : sv, sub_it
@OUTPUT     : tv
@RETURNS    : changed (1) or unchanged (0), based on whether voxels were thinned
@DESCRIPTION: Reduces one of the skeletonization subfield, based on the value 
              of sub_it
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */

int sub_reducer(Volume sv, Volume tv, int sub_it) {

  int changed = 0;
  int sizes[3];
  int v1, v2, v3;
  int shift;
  Real voxel_value;

  get_volume_sizes(sv, sizes);

  /* The algorithm does not thin voxels on the outer boundary of the volume,
     though these are used in the thinning decision procedure.  The if-else
     conditions below insure that depending on the subiteration the correct
     set of voxels are used in the thinning procedure.  (All voxels that 
     share a single diagonal edge.
  */
  for( v1 = 2;  v1 < sizes[0]-2;  v1++ ) {  
    for( v2 = 2;  v2 < sizes[1]-2;  v2++ ) {
      if (sub_it == 0) {
	if ((v1 % 2) == 0) {
	  if ((v2 % 2) == 0) 
	    shift = 0;
	  else
	    shift = 1;	  
	}
	else {
	  if ((v2 % 2) == 0) 
	    shift = 1;
	  else
	    shift = 0;	  
	}
      }
      else { //sub_it = 1
	if ((v1 % 2) == 0) {
	  if ((v2 % 2) == 0) 
	    shift = 1;
	  else
	    shift = 0;	  
	}
	else {
	  if ((v2 % 2) == 0) 
	    shift = 0;
	  else
	    shift = 1;	  
	}
      }
      for( v3 = 2;  v3 < sizes[2]-2;  v3++ ) {

	//shift insures right subfield is being thinned along each row.
	if (shift == (v3 % 2)) {	  
	  voxel_value = get_volume_real_value(sv, v1, v2, v3, 0, 0);

	  if (voxel_value > VALUE_0) {

	    set_volume_real_value(tv,v1,v2,v3,0,0,VALUE_1);

	    //check if voxel is end voxel, otherwise check whether it should
	    //be deleted
	    if (!end_voxel(sv,v1,v2,v3)){
	      if (deletable_not_preserving(sv,v1,v2,v3)) {	  
		set_volume_real_value(tv,v1,v2,v3,0,0,VALUE_0);
		changed = 1;
		continue;
	      }
	    }
	  }
	  else 
	    set_volume_real_value(tv,v1,v2,v3,0,0,VALUE_0);
	}	
	else {
	  voxel_value = get_volume_real_value(sv, v1, v2, v3, 0, 0);
	  set_volume_real_value(tv,v1,v2,v3,0,0,voxel_value);
	}
      }
    }
  }
  return changed;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : end_voxel
@INPUT      : vol, v1, v2, v3
@OUTPUT     : 
@RETURNS    : IS END VOXEL (1) or NOT END VOXEL (0)
@DESCRIPTION: Determines whether a voxel is or is not an end voxel.  (Non
              end voxels can be deleted.)  A voxel is an end voxel if ...
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */

int end_voxel(Volume vol,int v1, int v2, int v3) {
  
  int v1t,v2t,v3t;
  int adj = 0;
  Real voxel_value;


  for( v1t = -1;  v1t < 2;  v1t++ ) {  
    for( v2t = -1;  v2t < 2;  v2t++ ) {
      for( v3t = -1;  v3t < 2;  v3t++ ) {

	if (!((v1t == 0) && (v2t == 0) && (v3t == 0))) {	 
 	  voxel_value = get_volume_real_value(vol, v1+v1t, v2+v2t, v3+v3t, 0, 0);
	  if (voxel_value > VALUE_0) {
	    if (adj == 1)
	      return FALSE;
	    adj++;
	  }
	}

      } // v3t
    } // v2t
  } // v1t 
 
  return adj;
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : deletable_not_preserving
@INPUT      : sv, v1, v2, v3
@OUTPUT     : 
@RETURNS    : TRUE (1) or FALSE (0)
@DESCRIPTION: Determines whether a voxel has a neighbourhood which makes it
              deletable, but not one that also allows it to be preserved.
	      Please consult the reference given at the top of the file
	      for more clarification.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */

int deletable_not_preserving(Volume sv, int v1, int v2, int v3) {
  
  int deletable = FALSE;
  int order1[6] = { NONE, X90, X180, XN90, Y90, YN90 };
  int order2[12] = { NONE, Y90, Y180, YN90, XN90, XN90_ZN90,
		    XN90_Z90, XN90_Z180, Z90, Z90_X90,
                    Z90_X180, Z90_XN90};
  int i;

  //The conditions under which a voxel may vertically deletable  
  //must be checked under all possible rotations under which
  //the condition may hold.
  for (i = 0; i < 6; i++) {
    if (check_vd(sv,v1,v2,v3,order1[i])) {
      deletable = TRUE;
      if (check_vp(sv,v1,v2,v3,order1[i])) {
	deletable = FALSE;
	break;
      }
    }
  }

  if (deletable)
    return TRUE;

  //For the diagonally deletable case, there are 12 possible rotation
  //combinations which must be checked under which the condition may hold
  //Some  rotations are required to be width 3 deletable as well.
  for (i = 0; i < 12; i++) {
    if (check_dd(sv,v1,v2,v3,order2[i])) {
       if ((order2[i] == Z90_X90) ||
	  (order2[i] == Z90_X180) ||
	  (order2[i] == Y90) ||
	  (order2[i] == XN90_ZN90) ||
	  (order2[i] == XN90_Z180)) {
	if (check_dd_w3(sv,v1,v2,v3,order2[i]))
	  deletable = TRUE;
      }
      else
	deletable = TRUE;
      
      if (check_dp(sv,v1,v2,v3,order2[i])) {
	deletable = FALSE;
	break;
      }      
    }
  }

  return deletable;
} 
  

/* ----------------------------- MNI Header -----------------------------------
@NAME       : check_vd
@INPUT      : sv, v1, v2, v3, rotation
@OUTPUT     : 
@RETURNS    : TRUE (1) or FALSE (0)
@DESCRIPTION: Determines whether a voxel is vertically deletable, under the
              rotation specified.  The conditions for vd are specified
	      in the reference given above.	           
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */

int check_vd(Volume sv,int v1,int v2,int v3,rotation rot) {


  //the voxel above the voxel in question must be 0 and the voxel 
  //below must be 1
  if (!(!rotated_value(sv,v1,v2,v3,rot,1,0,0) &&
      rotated_value(sv,v1,v2,v3,rot,-1,0,0))) {
    return FALSE;
  }


  //under all rotations check whether in the midEW and midNS layers,
  //of the 26-adjacent cube to the current voxel, only one 1-component
  //is present.
  if (rotated_value(sv,v1,v2,v3,rot,1,-1,0)) {
    if (!rotated_value(sv,v1,v2,v3,rot,0,-1,0)) 
      return FALSE;
  }

  if (rotated_value(sv,v1,v2,v3,rot,1,1,0)) {
    if (!rotated_value(sv,v1,v2,v3,rot,0,1,0)) 
      return FALSE;
  }

  if (rotated_value(sv,v1,v2,v3,rot,1,0,1)) {
    if (!rotated_value(sv,v1,v2,v3,rot,0,0,1)) 
      return FALSE;
  }

  if (rotated_value(sv,v1,v2,v3,rot,1,0,-1)) {
    if (!rotated_value(sv,v1,v2,v3,rot,0,0,-1)) 
      return FALSE;
  }


  //This condition insures that if an upper corner of the 26-adjacent
  //cube is set, then of two other specific voxels is set.  (i.e., if
  //upper corner voxel (0,0,1) is set, then one of (0,1,0) or (1,0,0) 
  //must be set.
  if (rotated_value(sv,v1,v2,v3,rot,1,-1,1)) {
    if ((!rotated_value(sv,v1,v2,v3,rot,0,-1,0)) && 
	(!rotated_value(sv,v1,v2,v3,rot,0,0,1)))
      return FALSE;
  }  

  if (rotated_value(sv,v1,v2,v3,rot,1,1,1)) {
    if ((!rotated_value(sv,v1,v2,v3,rot,0,1,0)) && 
	(!rotated_value(sv,v1,v2,v3,rot,0,0,1)))
      return FALSE;
  }  

  if (rotated_value(sv,v1,v2,v3,rot,1,1,-1)) {
    if ((!rotated_value(sv,v1,v2,v3,rot,0,1,0)) && 
	(!rotated_value(sv,v1,v2,v3,rot,0,0,-1)))
      return FALSE;
  }  

  if (rotated_value(sv,v1,v2,v3,rot,1,-1,-1)) {
    if ((!rotated_value(sv,v1,v2,v3,rot,0,-1,0)) && 
	(!rotated_value(sv,v1,v2,v3,rot,0,0,-1)))
      return FALSE;
  }  

  return TRUE;

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : check_vp
@INPUT      : sv, v1, v2, v3, rotation
@OUTPUT     : 
@RETURNS    : TRUE (1) or FALSE (0)
@DESCRIPTION: If a voxel has been found to be deletable, conditions to preserve
              it may still be present.  This is the condition for the vertical
	      case.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */

int check_vp(Volume sv,int v1,int v2,int v3, rotation rot) {
  
  int adjacent[16][3];
  int traversed[24][3];
  int adj_it = 0;
  int at_voxel=0;
  int i,j,k;
  int tot_voxels=0;
  Real voxel_value;
  int n; 
  int found = 0;
  int absv1, absv2, absv3;
  int vmin1, vmin2, vmin3;
  int cnt = 0;
  int tk, tj, ti;

  int endinlayer = FALSE;
  
  //the current voxel must be an end voxel in the MID NS or MID EW layers.

  //check mid NS
  for (j = -1; j < 2; j++) {
    for (i = -1; i < 2 ; i++) {      
      if ((i==0) && (j==0))
	continue;
      if ((i==-1) && (j==0)) 
	continue;
      if (rotated_value(sv,v1,v2,v3,rot,j,0,i))
	endinlayer = TRUE;
    }
  }

  //check mid EW
  if (!endinlayer) {
    for (j = -1; j < 2; j++) {
      for (i = -1; i < 2 ; i++) {      
	if ((i==0) && (j==0))
	  continue;
	if ((i==-1) && (j==0)) 
	  continue;
	if (rotated_value(sv,v1,v2,v3,rot,0,j,i))
	  return FALSE;
      }
    }
  }
  
  //find all non-zero voxels adjacent to the current voxel, aside from the
  //ones whose states are already known, (the one directly above and
  //the one directly below)

  for (k = -1; k < 1; k++) {
    for (j = -1; j < 2 ; j++) {
      for (i = -1; i < 2 ; i++) {

	if ((j == 0) && (i == 0))
	  continue;
	  
	voxel_value = rotated_value(sv, v1, v2, v3, rot,k,j,i);

	if (voxel_value) {
	  adjacent[adj_it][0] = k;
	  adjacent[adj_it][1] = j;
	  adjacent[adj_it++][2] = i;
	}	    	
      } // i
    } // j
  } // k

  //if no other adjacent voxels are present then return false.
  if (!adj_it)
    return FALSE;

  traversed[0][0] = adjacent[0][0];
  traversed[0][1] = adjacent[0][1];
  traversed[0][2] = adjacent[0][2];
  tot_voxels = 1;
  at_voxel = 0;

  while (1) {
    for (k = -1; k < 2; k++) {
      for (j = -1; j < 2 ; j++) {
	for (i = -1; i < 2 ; i++) {
	  int valid = 1;

	  if ((k == 0) && (j == 0) && (i == 0))
	    continue;

	  if (((k+traversed[at_voxel][0]) > 1) || ((k+traversed[at_voxel][0]) < -1) ||
	      ((j+traversed[at_voxel][1]) > 1) || ((j+traversed[at_voxel][1]) < -1) ||
	      ((i+traversed[at_voxel][2]) > 1) || ((i+traversed[at_voxel][2]) < -1))
	    continue;

	  rotated_coords(sv,v1,v2,v3,rot,k+traversed[at_voxel][0],
			 j+traversed[at_voxel][1],i+traversed[at_voxel][2],
			 &absv1,&absv2,&absv3);
	  
	  rotated_coords(sv,v1, v2, v3,rot,-1,0,0,&vmin1,&vmin2,&vmin3);

	  if (((v1 == absv1)&&(v2 == absv2)&&(v3 == absv3))||
	      ((vmin1 == absv1)&&(vmin2 == absv2)&&(vmin3 == absv3)))
	    continue;

	  voxel_value = rotated_value(sv,v1,v2,v3,rot,k+traversed[at_voxel][0],
				      j+traversed[at_voxel][1],i+traversed[at_voxel][2]);

	  if (voxel_value == VALUE_0)
	    continue;


	  valid = 1;
	  for (n = 0; n < tot_voxels; n++) {
	    
	    if ((traversed[n][0] == (k+traversed[at_voxel][0])) &&
		(traversed[n][1] == (j+traversed[at_voxel][1])) &&
		(traversed[n][2] == (i+traversed[at_voxel][2]))) {
	      valid = 0;
	      break;
	    }	
	  } //n

	  if (!valid)
	    continue;

	  traversed[tot_voxels][0] = k+traversed[at_voxel][0];
	  traversed[tot_voxels][1] = j+traversed[at_voxel][1];
	  traversed[tot_voxels++][2] = i+traversed[at_voxel][2];

	} //i
      } //j
    } //k

    at_voxel++;
    if (at_voxel == tot_voxels) {
      break;
    }
  }

  for(i = 0; i < adj_it; i++) {
    found = 0;
    for(j = 0; j < tot_voxels; j++) {
      if ((adjacent[i][0] == traversed[j][0]) &&
	  (adjacent[i][1] == traversed[j][1]) &&
	  (adjacent[i][2] == traversed[j][2])) {
	found = 1;
	break;
      }
    }
    if (found == 0)
      return TRUE;
  }

  //if two adjacent voxels exist determine if the they, along with the
  //current and the lower neighbour voxel form a square.  If so the
  //voxel is preserved.
  if (adj_it == 2) {
    
    Real v_val1;
    Real v_val2;
    
    v_val1 = rotated_value(sv,v1,v2,v3,rot,0,-1,0);
    v_val2 = rotated_value(sv,v1,v2,v3,rot,-1,-1,0);
    
    if ((v_val1) && (v_val2))
      return TRUE;
    
    v_val1 = rotated_value(sv,v1,v2,v3,rot,0,1,0);
    v_val2 = rotated_value(sv,v1,v2,v3,rot,-1,1,0);
    
    if ((v_val1) && (v_val2))
      return TRUE;

    v_val1 = rotated_value(sv,v1,v2,v3,rot,0,0,1);
    v_val2 = rotated_value(sv,v1,v2,v3,rot,-1,0,1);
	
    if ((v_val1) && (v_val2))
      return TRUE;

    v_val1 = rotated_value(sv,v1,v2,v3,rot,0,0,-1);
    v_val2 = rotated_value(sv,v1,v2,v3,rot,-1,0,-1);
	
    if ((v_val1) && (v_val2))
      return TRUE;
  }

  return FALSE;
    
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : check_dd
@INPUT      : sv, v1, v2, v3, rotation
@OUTPUT     : 
@RETURNS    : TRUE (1) or FALSE (0)
@DESCRIPTION: The voxel in question is checked to see whether it is diagonally
              deletable. 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */

int check_dd(Volume sv,int v1,int v2,int v3,rotation rot) {

  int i,j;

  //the lw neighbour of the current voxel must be a 1.
  if (!rotated_value(sv,v1,v2,v3,rot,-1,0,-1)) {
    return FALSE;
  }


  //the current voxel must be an end voxel in the MID NS layer.
  for (j = -1; j < 2; j++) {
    for (i = -1; i < 2 ; i++) {      
      if ((i==0) && (j==0))
	continue;
      if ((i==-1) && (j==-1)) 
	continue;
      if (rotated_value(sv,v1,v2,v3,rot,j,0,i))
	return FALSE;
    }
  }

  //all voxels on the upper and east faces must be zero
  for (j = -1; j < 2; j++) {
    for (i = -1; i < 2 ; i++) {      
      if (rotated_value(sv,v1,v2,v3,rot,j,i,1))
	return FALSE;
      if (rotated_value(sv,v1,v2,v3,rot,1,j,i))
	return FALSE;
    }
  }

  return TRUE;

}

/*******************************************************************/

/* ----------------------------- MNI Header -----------------------------------
@NAME       : check_dp
@INPUT      : sv, v1, v2, v3, rotation
@OUTPUT     : 
@RETURNS    : TRUE (1) or FALSE (0)
@DESCRIPTION: The voxel in question is checked to see if it is diagonally 
              preserving.  
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */

int check_dp(Volume sv,int v1,int v2,int v3,rotation rot) {

  //The conditions for diagonal preservation are too complex to relate 
  //in ascii text.  They can be found in the reference given above.
  
  if (!(rotated_value(sv,v1,v2,v3,rot,-1,-1,0) || 
	rotated_value(sv,v1,v2,v3,rot,0,-1,-1) ||
	rotated_value(sv,v1,v2,v3,rot,-1,1,0) ||
	rotated_value(sv,v1,v2,v3,rot,0,1,-1))) {

    if (rotated_value(sv,v1,v2,v3,rot,0,-1,0) && 
	rotated_value(sv,v1,v2,v3,rot,-1,-1,-1))
      return TRUE;

    if (rotated_value(sv,v1,v2,v3,rot,0,1,0) && 
	rotated_value(sv,v1,v2,v3,rot,-1,1,-1))
      return TRUE;
  }
  
  
  if (!(rotated_value(sv,v1,v2,v3,rot,0,-1,0) || 
	rotated_value(sv,v1,v2,v3,rot,-1,-1,-1))) {
    if (!(rotated_value(sv,v1,v2,v3,rot,-1,-1,0) && 
	rotated_value(sv,v1,v2,v3,rot,0,-1,-1)))
      return TRUE;
  }
  
  if (!(rotated_value(sv,v1,v2,v3,rot,0,1,0) ||
	rotated_value(sv,v1,v2,v3,rot,-1,1,-1)))
    if (!(rotated_value(sv,v1,v2,v3,rot,-1,1,0) && 
	  rotated_value(sv,v1,v2,v3,rot,0,1,-1)))
      return TRUE;

  return FALSE;
}
  
/* ----------------------------- MNI Header -----------------------------------
@NAME       : check_dd_w3
@INPUT      : sv, v1, v2, v3, rotation
@OUTPUT     : 
@RETURNS    : TRUE (1) or FALSE (0)
@DESCRIPTION: The voxel in question is checked to determine whether it is 
              width 3 deletable.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */

int check_dd_w3(Volume sv,int v1,int v2,int v3,rotation rot) {
  
  //Refer to the paper mentioned above for information about this condition

  if (rotated_value(sv,v1,v2,v3,rot,0,0,-2) ||
      rotated_value(sv,v1,v2,v3,rot,-1,0,-2) ||
      rotated_value(sv,v1,v2,v3,rot,-2,0,-2) ||
      rotated_value(sv,v1,v2,v3,rot,-2,0,-1) ||
      rotated_value(sv,v1,v2,v3,rot,-2,0,0))
    return TRUE;
  
  return FALSE;  
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : reduce_ends
@INPUT      : sv
@OUTPUT     : tv
@RETURNS    : 
@DESCRIPTION: All end one voxels are removed.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */


void reduce_ends(Volume sv, Volume tv) {
  
  Real voxel_value;
  int num, v1, v2, v3;
  int zr, yr, xr;
  int sizes[3];
  get_volume_sizes(sv, sizes);
  
  //an end ones are deleted

  for( v1 = 1;  v1 < sizes[0]-1;  v1++ ) {  
    for( v2 = 1;  v2 < sizes[1]-1;  v2++ ) {  
      for( v3 = 1;  v3 < sizes[2]-1;  v3++ ) {  
	voxel_value = get_volume_real_value(sv, v1, v2, v3, 0, 0);
	if (voxel_value > VALUE_0) {
	  num = 0;
	  for( zr = -1;  zr < 2;  zr++ ) {  
	    for( yr = -1;  yr < 2;  yr++ ) {  
	      for( xr = -1;  xr < 2;  xr++ ) {  
		voxel_value = get_volume_real_value(sv, v1+zr, v2+yr, v3+xr, 0, 0);	    
		if (voxel_value > VALUE_0)
		  num++;
	      }
	    }
	  }  
	  if (num == 2) {
	    set_volume_real_value(tv,v1,v2,v3,0,0,VALUE_0);
	  }
	  else
	    set_volume_real_value(tv,v1,v2,v3,0,0,VALUE_1);
	}
	else
	    set_volume_real_value(tv,v1,v2,v3,0,0,VALUE_0);	  
      }
    }
  }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : rotated_value
@INPUT      : vol, v1, v2, v3, rot, zr, yr, xr
@OUTPUT     : 
@RETURNS    : value at voxel
@DESCRIPTION: Determine the value of the voxel at (v1,v2,v3) offset by (zr,yr,xr)
              with a rotation specified by rot around the v1,v2,v3 coordinate
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */
Real rotated_value(Volume vol, int v1, int v2, int v3, rotation rot,
		   int zr, int yr, int xr) {


  Real voxel_value;
  
  if (rot == NONE) {
    voxel_value = get_volume_real_value( vol, v1+zr, v2+yr, v3+xr, 0, 0);
    return voxel_value;
  }


  if (rot == X90) {
    voxel_value = get_volume_real_value( vol, v1+yr, v2-zr, v3+xr, 0, 0);
    return voxel_value;
  }
  if (rot == X180) {
    voxel_value = get_volume_real_value( vol, v1-zr, v2-yr, v3+xr, 0, 0);
    return voxel_value;
  }
  if (rot == XN90) {
    voxel_value = get_volume_real_value( vol, v1-yr, v2+zr, v3+xr, 0, 0);
    return voxel_value;
  }


  if (rot == Y90) {
    voxel_value = get_volume_real_value( vol, v1+xr, v2+yr, v3-zr, 0, 0);
    return voxel_value;
  }
  if (rot == Y180) {
    voxel_value = get_volume_real_value( vol, v1-zr, v2+yr, v3-xr, 0, 0);
    return voxel_value;
  }
  if (rot == YN90) {
    voxel_value = get_volume_real_value( vol, v1-xr, v2+yr, v3+zr, 0, 0);
    return voxel_value;
  }



  if (rot == Z90) {
    voxel_value = get_volume_real_value( vol, v1+zr, v2+xr, v3-yr, 0, 0);
    return voxel_value;
  }
  if (rot == Z180) {
    voxel_value = get_volume_real_value( vol, v1+zr, v2-yr, v3-xr, 0, 0);
    return voxel_value;
  }
   if (rot == ZN90) {
    voxel_value = get_volume_real_value( vol, v1+zr, v2-xr, v3+yr, 0, 0);
    return voxel_value;
  } 

  //DUAL TRANSFORMS

  if (rot == XN90_Z90) {
    voxel_value = get_volume_real_value( vol, v1-yr, v2+xr, v3-zr, 0, 0);
    return voxel_value;
  }
  if ((rot == XN90_ZN90)) {
    voxel_value = get_volume_real_value( vol, v1-yr, v2-xr, v3+zr, 0, 0);
    return voxel_value;
  }
  if ((rot == XN90_Z180)) {
    voxel_value = get_volume_real_value( vol, v1-yr, v2-zr, v3-xr, 0, 0);
    return voxel_value;
  }


  if ((rot == Z90_X90)) {
    voxel_value = get_volume_real_value( vol, v1+xr, v2-zr, v3-yr, 0, 0);
    return voxel_value;
  }
  if ((rot == Z90_XN90)) {
    voxel_value = get_volume_real_value( vol, v1-xr, v2+zr, v3-yr, 0, 0);
    return voxel_value;
  }
  if ((rot == Z90_X180)) {
    voxel_value = get_volume_real_value( vol, v1-zr, v2-xr, v3-yr, 0, 0);
    return voxel_value;
  }

  return VALUE_0;
  
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : rotated_value
@INPUT      : vol, v1, v2, v3, rot, zr, yr, xr
@OUTPUT     : absv1, absv2, absv3
@RETURNS    : 
@DESCRIPTION: Determine the absolute coordinates of the voxel after a rotation
              is applied centred at (v1, v2, v3) with offset (zr,yr,xr)
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */

void rotated_coords(Volume vol, int v1, int v2, int v3, rotation rot,
                   int zr, int yr, int xr, int* absv1, int* absv2, int* absv3) {
  

  if (rot == NONE) {
    *absv1 =  v1+zr; *absv2 =  v2+yr; *absv3 =  v3+xr;
    
  }


  if (rot == X90) {
    *absv1 =  v1+yr; *absv2 =  v2-zr; *absv3 =  v3+xr;
    
  }
  if (rot == X180) {
    *absv1 =  v1-zr; *absv2 =  v2-yr; *absv3 =  v3+xr;
    
  }
  if (rot == XN90) {
    *absv1 =  v1-yr; *absv2 =  v2+zr; *absv3 =  v3+xr;
    
  }


  if (rot == Y90) {
    *absv1 =  v1+xr; *absv2 =  v2+yr; *absv3 =  v3-zr;
    
  }
  if (rot == Y180) {
    *absv1 =  v1-zr; *absv2 =  v2+yr; *absv3 =  v3-xr;
    
  }
  if (rot == YN90) {
    *absv1 =  v1-xr; *absv2 =  v2+yr; *absv3 =  v3+zr;
    
  }



  if (rot == Z90) {
    *absv1 =  v1+zr; *absv2 =  v2+xr; *absv3 =  v3-yr;
    
  }
  if (rot == Z180) {
    *absv1 =  v1+zr; *absv2 =  v2-yr; *absv3 =  v3-xr;
    
  }
   if (rot == ZN90) {
    *absv1 =  v1+zr; *absv2 =  v2-xr; *absv3 =  v3+yr;
    
  } 

  //DUAL TRANSFORMS

  if (rot == XN90_Z90) {
    *absv1 =  v1-yr; *absv2 =  v2+xr; *absv3 =  v3-zr;
    
  }
  if ((rot == XN90_ZN90)) {
    *absv1 =  v1-yr; *absv2 =  v2-xr; *absv3 =  v3+zr;
    
  }
  if ((rot == XN90_Z180)) {
    *absv1 =  v1-yr; *absv2 =  v2-zr; *absv3 =  v3-xr;
    
  }


  if ((rot == Z90_X90)) {
    *absv1 =  v1+xr; *absv2 =  v2-zr; *absv3 =  v3-yr;
    
  }
  if ((rot == Z90_XN90)) {
    *absv1 =  v1-xr; *absv2 =  v2+zr; *absv3 =  v3-yr;
    
  }
  if ((rot == Z90_X180)) {
    *absv1 =  v1-zr; *absv2 =  v2-xr; *absv3 =  v3-yr;
    
  }


}
