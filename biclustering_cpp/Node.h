#ifndef CLUSTERING_NODE_H
#define CLUSTERING_NODE_H

#include <armadillo>
#include <map>
#include <set>
#include <vector>
#include <MatlabDataArray/TypedArray.hpp>

#define BRANCH_ON_U 0
#define BRANCH_ON_V 1

class Node {

public:

    std::vector<std::pair<int, int>> global_ml_pairs_U;
    std::vector<std::pair<int, int>> global_cl_pairs_U;
    std::vector<std::pair<int, int>> global_ml_pairs_V;
    std::vector<std::pair<int, int>> global_cl_pairs_V;

    // lower bound
    double lb;
    // upper bound
    double ub;
    // node id
    int id;

};


class SDPNode : public Node {


public:

    std::vector<arma::sp_mat> B_vector;

};

typedef struct NodeData {

    SDPNode *node;
    int i;
    int j;

} NodeData;

typedef struct JobData {

    int type;
    NodeData *node_data;

} JobData;

typedef struct SDPResult {

    int n; // current problem size
    int m; // current problem size
    double lb;
    double ub;
    arma::sp_mat Xu_assignment;
    arma::sp_mat Xv_assignment;
    int n_ineq;
    int cp_iter;
    int cp_flag;
    int branching_type;
    int branching_i;
    int branching_j;
    std::vector<arma::sp_mat> B_vector;

} SDPResult;


#endif //CLUSTERING_NODE_H
