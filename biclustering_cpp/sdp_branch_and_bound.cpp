#include <thread>
#include "matlab_util.h"
#include "sdp_branch_and_bound.h"
#include "JobQueue.h"
#include "util.h"
#include "config_params.h"
#include "Node.h"
#include "ThreadPool.h"

// root
SDPResult solve_sdp(std::unique_ptr<matlab::engine::MATLABEngine> &matlabPtr, matlab::data::ArrayFactory &factory, arma::mat &W, int k) {

    // convert data
    matlab::data::TypedArray<double> matlab_W = arma_to_matlab_matrix(factory, W);
    matlab::data::TypedArray<double> matlab_k = factory.createScalar<double>(k);

    // Create StructArray
    std::vector<std::string> f = {"n_threads", "bb_tol", "sdp_verbose", "sdp_tol", "cp_maxiter", "cp_tol", "cp_maxineq",
                                  "cp_percineq", "cp_epsineq", "cp_activeineq", "gurobi_verbose"};

    matlab::data::StructArray struct_matlab = factory.createStructArray({1, 1},f);
    struct_matlab[0]["n_threads"] = factory.createScalar<int>(matlab_session_threads_root);
    struct_matlab[0]["bb_tol"] = factory.createScalar<double>(branch_and_bound_tol);
    struct_matlab[0]["sdp_verbose"] = factory.createScalar<int>(sdp_solver_verbose);
    struct_matlab[0]["sdp_tol"] = factory.createScalar<double>(sdp_solver_tol);
    struct_matlab[0]["cp_maxiter"] = factory.createScalar<int>(cp_max_iter);
    struct_matlab[0]["cp_tol"] = factory.createScalar<double>(cp_tol);
    struct_matlab[0]["cp_maxineq"] = factory.createScalar<double>(cp_max_ineq);
    struct_matlab[0]["cp_percineq"] = factory.createScalar<double>(cp_perc_ineq);
    struct_matlab[0]["cp_epsineq"] = factory.createScalar<double>(cp_eps_ineq);
    struct_matlab[0]["cp_activeineq"] = factory.createScalar<double>(cp_eps_active);
    struct_matlab[0]["gurobi_verbose"] = factory.createScalar<double>(heuristic_verbose);

    std::vector<matlab::data::Array> args({matlab_W, matlab_k, struct_matlab});

    // Call MATLAB function and return result
    const size_t n_return = 1;
    matlabPtr->eval(u"clear");
    std::vector<matlab::data::Array> result = matlabPtr->feval(u"call_solve_biclustering_root",n_return, args);

    matlab::data::StructArray structArray = result[0];
    matlab::data::TypedArray<double> field_best_lb = structArray[0]["best_lb"];
    matlab::data::TypedArray<double> field_best_ub = structArray[0]["best_ub"];
    matlab::data::SparseArray<double> field_best_Xu = structArray[0]["best_Xu"];
    matlab::data::SparseArray<double> field_best_Xv = structArray[0]["best_Xv"];
    matlab::data::TypedArray<double> field_cp_iter = structArray[0]["cp_iter"];
    matlab::data::TypedArray<double> field_cp_flag = structArray[0]["cp_flag"];
    matlab::data::CellArray field_best_B_cell = structArray[0]["best_B_cell"];
    matlab::data::TypedArray<double> field_i_idx = structArray[0]["i_idx"];
    matlab::data::TypedArray<double> field_j_idx = structArray[0]["j_idx"];
    matlab::data::TypedArray<double> field_branching_type = structArray[0]["branching_type"];

    double lb = field_best_lb[0];
    double ub = field_best_ub[0];
    arma::sp_mat Xu_assignment = matlab_to_arma_sparse(field_best_Xu);
    arma::sp_mat Xv_assignment = matlab_to_arma_sparse(field_best_Xv);
    int cp_iter = (int) field_cp_iter[0];
    int cp_flag = (int) field_cp_flag[0];
    std::vector<arma::sp_mat> B_vector = matlab_to_arma_vector_sp_mat(field_best_B_cell);
    int n_ineq = (int) B_vector.size();
    int i_idx = (int) field_i_idx[0];
    int j_idx = (int) field_j_idx[0];
    int branching_type = (int) field_branching_type[0];
    int n = (int) W.n_rows;
    int m = (int) W.n_cols;

    return SDPResult{n, m, lb, ub, Xu_assignment, Xv_assignment, n_ineq, cp_iter, cp_flag, branching_type, i_idx, j_idx, B_vector};
}

// lower bound must link and cannot link
SDPResult solve_sdp(std::unique_ptr<matlab::engine::MATLABEngine> &matlabPtr, matlab::data::ArrayFactory &factory,
        arma::mat &W, int k, std::vector<arma::sp_mat> &parent_B_vector, double global_lb, arma::sp_mat &global_Xu, arma::sp_mat &global_Xv,
        std::vector<std::pair<int, int>> &global_ml_pairs_U, std::vector<std::pair<int, int>> &global_cl_pairs_U,
        std::vector<std::pair<int, int>> &global_ml_pairs_V, std::vector<std::pair<int, int>> &global_cl_pairs_V) {

    // convert data
    matlab::data::TypedArray<double> matlab_W = arma_to_matlab_matrix(factory, W);
    matlab::data::TypedArray<double> matlab_k = factory.createScalar<double>(k);
    matlab::data::TypedArray<double> matlab_ml_U = vector_pair_to_matlab_matrix(factory, global_ml_pairs_U);
    matlab::data::TypedArray<double> matlab_cl_U = vector_pair_to_matlab_matrix(factory, global_cl_pairs_U);
    matlab::data::TypedArray<double> matlab_ml_V = vector_pair_to_matlab_matrix(factory, global_ml_pairs_V);
    matlab::data::TypedArray<double> matlab_cl_V = vector_pair_to_matlab_matrix(factory, global_cl_pairs_V);

    matlab::data::CellArray matlab_B_cell = arma_to_matlab_cell(factory, parent_B_vector);
    matlab::data::SparseArray<double> matlab_global_Xu = arma_to_matlab_sparse(factory, global_Xu);
    matlab::data::SparseArray<double> matlab_global_Xv = arma_to_matlab_sparse(factory, global_Xv);

    matlab::data::TypedArray<double> matlab_glb = factory.createScalar<double>(global_lb);

    // Create StructArray
    std::vector<std::string> f = {"n_threads", "bb_tol", "sdp_verbose", "sdp_tol", "cp_maxiter", "cp_tol",
                                  "cp_maxineq", "cp_percineq", "cp_epsineq", "cp_activeineq", "gurobi_verbose"};

    matlab::data::StructArray struct_matlab = factory.createStructArray({1, 1},f);
    struct_matlab[0]["n_threads"] = factory.createScalar<int>(matlab_session_threads_child);
    struct_matlab[0]["bb_tol"] = factory.createScalar<double>(branch_and_bound_tol);
    struct_matlab[0]["sdp_verbose"] = factory.createScalar<int>(sdp_solver_verbose);
    struct_matlab[0]["sdp_tol"] = factory.createScalar<double>(sdp_solver_tol);
    struct_matlab[0]["cp_maxiter"] = factory.createScalar<int>(cp_max_iter);
    struct_matlab[0]["cp_tol"] = factory.createScalar<double>(cp_tol);
    struct_matlab[0]["cp_maxineq"] = factory.createScalar<double>(cp_max_ineq);
    struct_matlab[0]["cp_percineq"] = factory.createScalar<double>(cp_perc_ineq);
    struct_matlab[0]["cp_epsineq"] = factory.createScalar<double>(cp_eps_ineq);
    struct_matlab[0]["cp_activeineq"] = factory.createScalar<double>(cp_eps_active);
    struct_matlab[0]["gurobi_verbose"] = factory.createScalar<double>(heuristic_verbose);

    std::vector<matlab::data::Array> args({matlab_W, matlab_k, matlab_ml_U, matlab_cl_U, matlab_ml_V, matlab_cl_V,
                                           matlab_B_cell, matlab_glb, matlab_global_Xu, matlab_global_Xv, struct_matlab});

    // Call MATLAB function and return result
    const size_t n_return = 1;
    matlabPtr->eval(u"clear");
    std::vector<matlab::data::Array> result = matlabPtr->feval(u"call_solve_biclustering_child",n_return, args);

    matlab::data::StructArray structArray = result[0];
    matlab::data::TypedArray<double> field_best_lb = structArray[0]["best_lb"];
    matlab::data::TypedArray<double> field_best_ub = structArray[0]["best_ub"];
    matlab::data::SparseArray<double> field_best_Xu = structArray[0]["best_Xu"];
    matlab::data::SparseArray<double> field_best_Xv = structArray[0]["best_Xv"];
    matlab::data::TypedArray<double> field_cp_iter = structArray[0]["cp_iter"];
    matlab::data::TypedArray<double> field_cp_flag = structArray[0]["cp_flag"];
    matlab::data::CellArray field_best_B_cell = structArray[0]["best_B_cell"];
    matlab::data::TypedArray<double> field_n = structArray[0]["n"];
    matlab::data::TypedArray<double> field_m = structArray[0]["m"];
    matlab::data::TypedArray<double> field_i_idx = structArray[0]["i_idx"];
    matlab::data::TypedArray<double> field_j_idx = structArray[0]["j_idx"];
    matlab::data::TypedArray<double> field_branching_type = structArray[0]["branching_type"];

    double lb = field_best_lb[0];
    double ub = field_best_ub[0];
    arma::sp_mat Xu_assignment = matlab_to_arma_sparse(field_best_Xu);
    arma::sp_mat Xv_assignment = matlab_to_arma_sparse(field_best_Xv);
    int cp_iter = (int) field_cp_iter[0];
    int cp_flag = (int) field_cp_flag[0];
    std::vector<arma::sp_mat> B_vector = matlab_to_arma_vector_sp_mat(field_best_B_cell);
    int n_ineq = (int) B_vector.size();
    int branching_type = (int) field_branching_type[0];
    int i_idx = (int) field_i_idx[0];
    int j_idx = (int) field_j_idx[0];
    int n = (int) field_n[0];
    int m = (int) field_m[0];

    return SDPResult{n, m, lb, ub, Xu_assignment, Xv_assignment, n_ineq, cp_iter, cp_flag, branching_type, i_idx, j_idx, B_vector};
}



std::pair<JobData *, JobData *> create_cl_ml_jobs(double node_gap, SDPNode *node, int branching_type, int i, int j, NodeData *parent, SharedData *shared_data) {

	if (node_gap <= branch_and_bound_tol) {
        delete(node);
        if (parent != nullptr) {
            delete (parent->node);
            delete (parent);
        }
        return std::make_pair(nullptr, nullptr);
    }

    if (i == -1 && j == -1) {

        const std::lock_guard<std::mutex> lock(shared_data->queueMutex);

        log_file << "PRUNING BY OPTIMALITY " << node->id << "\n";
        delete (node);
        if (parent != nullptr) {
            delete (parent->node);
            delete (parent);
        }
        return std::make_pair(nullptr, nullptr);

        // mutex is automatically released when lock goes out of scope
    }

    auto *cl_data = new NodeData();
    cl_data->node = new SDPNode(*node);
    cl_data->i = i;
    cl_data->j = j;

    auto *ml_data = new NodeData();
    ml_data->node = new SDPNode(*node);
    ml_data->i = i;
    ml_data->j = j;

    auto *ml_job_data = new JobData();
    auto *cl_job_data = new JobData();

    switch (branching_type) {
        case BRANCH_ON_U:
            ml_job_data->type = MUST_LINK_U;
            ml_job_data->node_data = ml_data;
            cl_job_data->type = CANNOT_LINK_U;
            cl_job_data->node_data = cl_data;
            break;
        case BRANCH_ON_V:
            ml_job_data->type = MUST_LINK_V;
            ml_job_data->node_data = ml_data;
            cl_job_data->type = CANNOT_LINK_V;
            cl_job_data->node_data = cl_data;
            break;
        default:
            std::cout << "Unknown branching type" << "\n";
            exit(EXIT_FAILURE);
    }

    if (parent != nullptr) {
        delete (parent->node);
        delete (parent);
    }

    delete (node);

    return std::make_pair(cl_job_data, ml_job_data);

}

std::pair<JobData *, JobData *> build_child_problem(int job_type, NodeData *node_data, InputData *input_data, SharedData  *shared_data) {

    auto *matlab_struct = new MatlabStruct();
    matlab_struct->matlabPtr = start_matlab(sdp_solver_folder, gurobi_folder);

    // generate cannot link child
    auto child_node = new SDPNode();
    auto parent = node_data->node;

    double parent_gap = (parent->ub - shared_data->global_lb) / shared_data->global_lb;
    if (parent_gap <= branch_and_bound_tol) return std::make_pair(nullptr, nullptr);

    child_node->global_ml_pairs_U = parent->global_ml_pairs_U;
    child_node->global_ml_pairs_V = parent->global_ml_pairs_V;
    child_node->global_cl_pairs_U = parent->global_cl_pairs_U;
    child_node->global_cl_pairs_V = parent->global_cl_pairs_V;

    switch (job_type) {
        case MUST_LINK_U:
            child_node->global_ml_pairs_U.emplace_back(node_data->i, node_data->j);
            break;
        case CANNOT_LINK_U:
            child_node->global_cl_pairs_U.emplace_back(node_data->i, node_data->j);
            break;
        case MUST_LINK_V:
            child_node->global_ml_pairs_V.emplace_back(node_data->i, node_data->j);
            break;
        case CANNOT_LINK_V:
            child_node->global_cl_pairs_V.emplace_back(node_data->i, node_data->j);
            break;
        default:
            std::cout << "Unknown job type" << "\n";
            exit(EXIT_FAILURE);
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    SDPResult  sdp_result = solve_sdp(matlab_struct->matlabPtr, matlab_struct->factory,
                                      input_data->W, input_data->k, parent->B_vector,
                                      shared_data->global_lb, shared_data->global_Xu, shared_data->global_Xv,
                                      child_node->global_ml_pairs_U, child_node->global_cl_pairs_U,
                                      child_node->global_ml_pairs_V, child_node->global_cl_pairs_V);

    int n = sdp_result.n;
    int m = sdp_result.m;
    int cp_iter = sdp_result.cp_iter;
    int cp_flag = sdp_result.cp_flag;
    int n_ineq = sdp_result.n_ineq;
    child_node->lb = std::max(sdp_result.lb, parent->lb);
    child_node->ub = sdp_result.ub;
    child_node->B_vector = sdp_result.B_vector;
    arma::sp_mat Xu_assignment = sdp_result.Xu_assignment;
    arma::sp_mat Xv_assignment = sdp_result.Xv_assignment;
    int branching_type = sdp_result.branching_type;
    int i = sdp_result.branching_i;
    int j = sdp_result.branching_j;

    double node_gap;

    {
        const std::lock_guard<std::mutex> lock(shared_data->queueMutex);

        bool ub_updated = false;
        if (shared_data->global_lb - child_node->lb <= -branch_and_bound_tol) {
            // update global lower bound
            shared_data->global_lb = child_node->lb;
            shared_data->global_Xu = Xu_assignment;
            shared_data->global_Xv = Xv_assignment;
            ub_updated = true;
        }

        child_node->id = shared_data->n_nodes;
        shared_data->n_nodes++;
        int open = shared_data->queue->getSize();

        node_gap = (child_node->ub - shared_data->global_lb) / child_node->ub;
        double gap = node_gap;
        Node *max_ub_node = shared_data->queue->getMaxUb();
        if (max_ub_node != nullptr)
            gap = (max_ub_node->ub - shared_data->global_lb) / max_ub_node->ub;

        shared_data->gap = gap;

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        auto time = (double) duration.count();

        print_log_sdp(n, m,parent->id, child_node->id, parent->ub, child_node->ub, time, cp_iter, cp_flag, n_ineq,
                      child_node->lb, shared_data->global_lb,branching_type, node_data->i, node_data->j, node_gap, shared_data->gap, open, ub_updated);

    }

    delete(matlab_struct);

    return create_cl_ml_jobs(node_gap, child_node, branching_type, i, j, node_data, shared_data);
}


std::pair<JobData *, JobData *> build_root_problem(MatlabStruct *matlab_struct, InputData *input_data, SharedData *shared_data) {

    // init root
    SDPNode *root;
    root = new SDPNode();
    root->id = shared_data->n_nodes;
    root->global_ml_pairs_U = {};
    root->global_cl_pairs_U = {};
    root->global_ml_pairs_V = {};
    root->global_cl_pairs_V = {};

    auto start_time = std::chrono::high_resolution_clock::now();

    SDPResult sdp_result = solve_sdp(matlab_struct->matlabPtr, matlab_struct->factory, input_data->W, input_data->k);

    int n = sdp_result.n;
    int m = sdp_result.m;
    int cp_iter = sdp_result.cp_iter;
    int cp_flag = sdp_result.cp_flag;
    int n_ineq = sdp_result.n_ineq;
    root->lb = sdp_result.lb;
    root->ub = sdp_result.ub;
    root->B_vector = sdp_result.B_vector;
    arma::sp_mat Xu_assignment = sdp_result.Xu_assignment;
    arma::sp_mat Xv_assignment = sdp_result.Xv_assignment;
    int branching_type = sdp_result.branching_type;
    int i = sdp_result.branching_i;
    int j = sdp_result.branching_j;

    // update shared data
    shared_data->global_lb = root->lb;
    shared_data->global_Xu = Xu_assignment;
    shared_data->global_Xv = Xv_assignment;
    shared_data->n_nodes++;

    int open = shared_data->queue->getSize();

    double node_gap = (root->ub - shared_data->global_lb) / root->ub;
    shared_data->gap = node_gap;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    auto time = (double) duration.count();

    print_log_sdp(n, m,-1, root->id, std::numeric_limits<double>::infinity(),
                  root->ub, time, cp_iter, cp_flag, n_ineq, root->lb,shared_data->global_lb,
                  branching_type, -1, -1, node_gap, node_gap, open, true);

    //return std::make_pair(nullptr, nullptr);
    return create_cl_ml_jobs(node_gap, root, branching_type, i, j, nullptr, shared_data);
}

void sdp_branch_and_bound(arma::mat &W, int k) {

    int n_thread = branch_and_bound_parallel;

    JobAbstractQueue *queue;
    switch (branch_and_bound_visiting_strategy) {
        case DEPTH_FIRST:
            queue = new JobStack();
            break;
        case BEST_FIRST:
            queue = new JobPriorityQueue();
            break;
        case BREADTH_FIRST:
            queue = new JobQueue();
            break;
        default:
            queue = new JobPriorityQueue();
    }

    auto *shared_data = new SharedData();
    shared_data->global_lb = -std::numeric_limits<double>::infinity();
    shared_data->n_nodes = 0;
    shared_data->queue = queue;

    shared_data->threadStates.reserve(n_thread);
    for (int i = 0; i < n_thread; i++) {
        shared_data->threadStates.push_back(false);
    }

    auto *input_data = new InputData();
    input_data->W = W; // data matrix
    input_data->k = k; // number of clusters

    ThreadPool pool(shared_data, input_data, n_thread);
    
    print_header_sdp();

    auto start_all = std::chrono::high_resolution_clock::now();
    
    auto *matlab_struct = new MatlabStruct();
    matlab_struct->matlabPtr = start_matlab(sdp_solver_folder, gurobi_folder);

    std::pair<JobData *, JobData *> jobs = build_root_problem(matlab_struct, input_data, shared_data);

    delete (matlab_struct);
    
    double root_gap = shared_data->gap;

    JobData *cl_job = jobs.first;
    JobData *ml_job = jobs.second;
    if (cl_job != nullptr && ml_job != nullptr) {
        pool.addJob(cl_job);
        pool.addJob(ml_job);
    }

    while (true) {

        {
            std::unique_lock<std::mutex> l(shared_data->queueMutex);
            while (is_thread_pool_working(shared_data->threadStates) && shared_data->n_nodes < branch_and_bound_max_nodes) {
                shared_data->mainConditionVariable.wait(l);
            }

            if (shared_data->queue->empty() || shared_data->n_nodes >= branch_and_bound_max_nodes)
                break;
        }

    }

    auto end_all = std::chrono::high_resolution_clock::now();
    auto duration_all = std::chrono::duration_cast<std::chrono::seconds>(end_all - start_all);

    pool.quitPool();

    if (queue->empty()) shared_data->gap = 0.0;

    log_file << "\n";
    log_file << "TIME: " << duration_all.count() << " sec\n";
    log_file << "NODES: " << shared_data->n_nodes << "\n";
    log_file << "ROOT_GAP: " << std::max(0.0, root_gap) << "\n";
    log_file << "GAP: " << std::max(0.0, shared_data->gap) << "\n";
    log_file << "OPT: " << shared_data->global_lb << "\n\n";

    // normalized partition matrices
    arma::sp_mat result_Xu = shared_data->global_Xu;
    arma::sp_mat result_Xv = shared_data->global_Xv;
    save_X_to_file(result_Xu, result_Xv);

    // free memory
    delete (input_data);
    delete (queue);
    delete (shared_data);

}
