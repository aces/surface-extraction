#ifndef  DEF_SUBDIV_PROTOTYPES
#define  DEF_SUBDIV_PROTOTYPES

public  void  initialize_uniform_subdiv(
    uniform_subdiv_struct  *subdiv,
    Point                  *min_limits,
    Point                  *max_limits,
    int                    nx,
    int                    ny,
    int                    nz );

public  void  delete_uniform_subdiv(
    uniform_subdiv_struct  *subdiv );

public  void  add_point_to_uniform_subdiv(
    uniform_subdiv_struct  *subdiv,
    Point                  *point,
    int                    point_index );

public  void  remove_point_from_uniform_subdiv(
    uniform_subdiv_struct  *subdiv,
    Point                  *point,
    int                    point_index );

public  void  move_point_in_uniform_subdiv(
    uniform_subdiv_struct  *subdiv,
    Point                  *prev_point,
    Point                  *new_point,
    int                    point_index );

public  void  create_uniform_subdiv(
    uniform_subdiv_struct  *subdiv,
    int                    n_points,
    Point                  points[],
    int                    n_cells );

public  int  get_subdiv_points_near_point(
    uniform_subdiv_struct  *subdiv,
    Point                  *point,
    Real                   distance,
    int                    points_list[] );
#endif
