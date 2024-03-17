function result = call_solve_biclustering_child(W, k, init_ML_U, init_CL_U, init_ML_V, init_CL_V, init_B_cell, global_lb, global_Xu, global_Xv, params)
    
    disp(init_ML_U)
    disp(init_CL_U)
    disp(init_ML_V)
    disp(init_CL_V)

    [n, m] = size(W);
    result = solve_biclustering_shrinking(W, k, init_ML_U, init_CL_U, init_ML_V, init_CL_V, init_B_cell, global_lb, global_Xu, global_Xv, params);
    result.best_lb = max(result.lb_list);
    if result.cp_flag == 1
        result.branching_type = -1;
        result.i_idx = -1;
        result.j_idx = -1;
    else
        [branching_type, i_idx, j_idx, ~] = get_branching_pair(result.best_Z, n, m);
        result.branching_type = branching_type;
        result.i_idx = i_idx;
        result.j_idx = j_idx;
    end
    
    %disp(result)
    
end