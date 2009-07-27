/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef  DEF_SUBDIV_H
#define  DEF_SUBDIV_H

#include  <bicpl.h>

typedef  struct
{
    int    n_points;
    int    *points;
}
subdiv_cell;

typedef  struct
{
    Point         limits[2];
    int           sizes[N_DIMENSIONS];
    subdiv_cell   ***cells;
}
uniform_subdiv_struct;

#include  <subdiv_prototypes.h>

#endif
