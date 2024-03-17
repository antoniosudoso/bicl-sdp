#ifndef CLUSTERING_UTIL_H
#define CLUSTERING_UTIL_H

#include <map>
#include <set>
#include <vector>
#include <armadillo>

void save_X_to_file(arma::sp_mat &Xu, arma::sp_mat &Xv);
/*
void print_header_sdp(std::ostream &log_file);
void print_log_sdp(std::ostream &log_file, int n, int m, int node_parent, int node, double ub_parent, double ub,
                   double time, int cp_iter, int cp_flag, int n_ineq, double lb, double glb,
                   int type, int i, int j, double node_gap, double gap, int open, bool update);
*/
void print_header_sdp();
void print_log_sdp(int n, int m, int node_parent, int node, double ub_parent, double ub,
                   double time, int cp_iter, int cp_flag, int n_ineq, double lb, double glb,
                   int type, int i, int j, double node_gap, double gap, int open, bool update);
#endif //CLUSTERING_UTIL_H
