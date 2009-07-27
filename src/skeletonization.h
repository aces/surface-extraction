/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/


//maximum number of subiterations to do for each iteration of algorithm
#define MAX_SUB 2

#define VALUE_0 0
#define VALUE_1 1

//typedefs used to rotate voxel coordinate system to check symmetric
//conditions
typedef enum { NONE, X90, X180, XN90, Y90, Y180, YN90,
	       Z90, Z180, ZN90, XN90_Z90, XN90_ZN90, XN90_Z180,
	       Z90_X90, Z90_XN90, Z90_X180 } rotation;

int sub_reducer(Volume vol1, Volume vol2, int sub_it);

int end_voxel(Volume sv, int v1, int v2, int v3);

int deletable_not_preserving(Volume sv, int v1, int v2, int v3);

int check_vd(Volume sv,int v1,int v2,int v3,rotation rot);

int check_vp(Volume sv, int v1, int v2, int v3,rotation rot);

int check_dd(Volume sv,int v1,int v2,int v3,rotation rot);

int check_dp(Volume sv,int v1,int v2, int v3,rotation rot);

int check_dd_w3(Volume sv,int v1,int v2,int v3,rotation rot);

void reduce_ends(Volume sv, Volume tv);

Real rotated_value(Volume vol, int v1, int v2, int v3, rotation rot,
		   int xr, int yr, int zr);

void reducer(Volume vol1, Volume vol2);

void rotated_coords(Volume vol, int v1, int v2, int v3, rotation rot,
                   int zr, int yr, int xr, int *absv1, int *absv2, int *absv3);
