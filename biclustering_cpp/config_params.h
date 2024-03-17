#ifndef CLUSTERING_CONFIG_PARAMS_H
#define CLUSTERING_CONFIG_PARAMS_H

#define MUST_LINK_U 0
#define CANNOT_LINK_U 1
#define MUST_LINK_V 2
#define CANNOT_LINK_V 3

#define BEST_FIRST 0
#define DEPTH_FIRST 1
#define BREADTH_FIRST 2

// data full path
extern const char *data_path;
extern const char *log_path;
extern const char *result_path;
extern std::ofstream log_file;

// branch and bound
extern double branch_and_bound_tol;
extern int branch_and_bound_parallel;
extern int branch_and_bound_max_nodes;
extern int branch_and_bound_visiting_strategy;

// matlab
extern int matlab_session_threads_root;
extern int matlab_session_threads_child;

// sdpnal
extern const char *sdp_solver_folder;
extern double sdp_solver_tol;
extern int sdp_solver_verbose;

// cutting plane
extern int cp_max_iter;
extern double cp_tol;
extern int cp_max_ineq;
extern double cp_perc_ineq;
extern double cp_eps_ineq;
extern double cp_eps_active;

// heuristic
extern const char *gurobi_folder;
extern int heuristic_verbose;

#endif //CLUSTERING_CONFIG_PARAMS_H
