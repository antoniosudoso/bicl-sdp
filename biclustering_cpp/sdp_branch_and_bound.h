#ifndef CLUSTERING_SDP_BRANCH_AND_BOUND_H
#define CLUSTERING_SDP_BRANCH_AND_BOUND_H

#include <armadillo>
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"
#include "JobQueue.h"

typedef struct MatlabStruct {

    std::unique_ptr<matlab::engine::MATLABEngine> matlabPtr;
    matlab::data::ArrayFactory factory;

} MatlabStruct;


typedef struct SharedData {

    // Between workers and main
    std::condition_variable mainConditionVariable;
    std::vector<bool> threadStates;

    // Queue of requests waiting to be processed
    JobAbstractQueue *queue;
    // This condition variable is used for the threads to wait until there is work to do
    std::condition_variable queueConditionVariable;
    // Mutex to protect queue
    std::mutex queueMutex;

    double global_lb;
    arma::sp_mat global_Xu;
    arma::sp_mat global_Xv;
    double gap;
    int n_nodes;

} SharedData;

typedef struct InputData {

    arma::mat W;
    int k;

} InputData;


void sdp_branch_and_bound(arma::mat &W, int k);
std::pair<JobData *, JobData *> build_root_problem(MatlabStruct *matlab_struct, InputData *input_data, SharedData *shared_data);
std::pair<JobData *, JobData *> build_child_problem(int job_type, NodeData *job_data, InputData *input_data, SharedData *shared_data);

#endif //CLUSTERING_SDP_BRANCH_AND_BOUND_H
