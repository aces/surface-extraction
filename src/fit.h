#ifndef  _DEF_FIT_H
#define  _DEF_FIT_H

#include  <volume_io/internal_volume_io.h>
#include  <deform.h>
#include  <subdiv.h>
#include  <special_geometry.h>

#define  INVALID_ID     -1

#define  INVALID_CONSTRAINT_VALUE    1.0e10

#define  NO_MAX_VALUE   1.0e30

typedef  struct
{
    object_struct           *object;
    int                     n_points;
    Point                   *points;
    Vector                  plane_normal;
    BOOLEAN                 subdiv_present;

    int                     *n_neighbours;
    int                     **neighbours;
    int                     total_neighbours;
    int                     *total_neighbours_list;
}
surface_struct;

typedef   struct
{
    int                       *n_neighbours;
    int                       **neighbours;
    int                       total_neighbours;
    int                       *total_neighbours_list;
    Real                      *total_lengths;
    Real                      *lengths;
    Real                      *curvatures;
    Real                      stretch_scale;
    Real                      min_stretch_factor, max_stretch_factor;
    Real                      stretch_weight;
    Real                      min_curvature_offset, max_curvature_offset;
    Real                      curvature_weight;
} model_info_struct;

typedef  enum  { SELF_INTERSECT, SURFACE_SURFACE,
                 STATIC_POINTS, DYNAMIC_POINTS }
               Constraint_types;

typedef  struct
{
    int       surface_index;
    Real      min_distance1;
    Real      min_distance2;
    Real      weight;
} self_intersect_struct;

typedef  struct
{
    int       surface_index1;
    int       surface_index2;
    Real      min_distance1;
    Real      min_distance2;
    Real      weight;
} surface_surface_struct;

typedef  struct
{
    int            surface_index1;
    int            surface_index2;
    int            *indices1;
    int            *indices2;
    float          *desired_distances;
    Real           min_distance_ratio;
    Real           max_distance_ratio;
    Real           weight;
} dynamic_points_struct;

typedef  struct
{
    int            surface_index;
    Smallest_int   *constrained_flags;
    Real           min_distance;
    Real           desired_distance;
    Real           max_distance;
    Real           weight;
    Point          *points;
} static_points_struct;

typedef   struct
{
    Constraint_types  type;
    union
    {
    self_intersect_struct    self_intersect;
    surface_surface_struct   surface_surface;
    dynamic_points_struct    dynamic_points;
    static_points_struct     static_points;
    } specific;
} constraint_struct;

typedef  enum  { THRESHOLD_VALUE, SURFACE_DISTANCE } Surface_types;

typedef struct
{
    Real    weight;
    Real    diameter;
    Volume  volume;
    BOOLEAN use_derivs[3];
    float   **features;
} multi_scale_struct;

typedef  struct
{
    Volume                       volume;
    Volume                       label_volume;
    voxel_coef_struct            voxel_lookup;
    bitlist_3d_struct            done_bits;
    bitlist_3d_struct            surface_bits;
    Surface_types                surface_type;
    boundary_definition_struct   boundary;
    Real                         max_outwards_distance;
    Real                         max_inwards_distance;
    Real                         weight;
    BOOLEAN                      threshold_per_vertex;
    Real                         *thresholds;
    Real                         max_diff;
    Normal_directions            normal_direction;

    Real                         multi_scale_weight;
    int                          n_multi_scales;
    multi_scale_struct           *scales;

    Real                         oversample_weight;
    int                          n_oversamples;
} surface_data_struct;

typedef  struct
{
    int                 n_surfaces;
    surface_struct      *surfaces;
    int                 *n_models;
    model_info_struct   **models;
    surface_data_struct *data;
    int                 n_constraints;
    constraint_struct   *constraints;
}
Deform_struct;

typedef struct
{
    Real      value;
    Real      deriv_x;
    Real      deriv_y;
    Real      deriv_z;
    BOOLEAN   found;
    Point     boundary;
} edge_info_struct;

#define  MAX_MODELS   3

typedef struct
{
    Real   volume_term;
    int    n_models;
    Real   stretch_term[MAX_MODELS];
    Real   curvature_term[MAX_MODELS];
    Real   constraint_term;
} fit_report_struct;

typedef  enum  {  GRID_OF_POINTS, FAST_GRID_OF_POINTS,
                  SPHERE_OF_POINTS } Grid_types;

#include  <fit_prototypes.h>

#endif
