function [result] = solve_biclustering_shrinking(W, k, init_ML_U, init_CL_U, init_ML_V, init_CL_V, init_B_cell, global_lb, global_Xu, global_Xv, params)

    % SOLVE SDP RELAXATION FOR BICLUSTERING
    
    % W: weight matrix (n x m)
    % k: number of clusters
    % init_ML_U: matrix of must-link constraints between vertices in U (n_ml_U x 2)
    % init_CL_U: matrix of cannot-link constraints between vertices in U (n_cl_U x 2)
    % init_ML_V: matrix of must-link constraints between vertices in V (n_ml_V x 2)
    % init_CL_V: matrix of cannot-link constraints between vertices in V (n_cl_V x 2)
    % init_B_cell: cell array of initial inequalities (for empty pass
    % cell(0)). Note that, inequalities are stored with original indices
    % global_lb: best known lower bound (double)
    % global_Xu: best known assignment 0-1 matrix (n x k)
    % global_Xv: best known assignment 0-1 matrix (m x k)
    
    % params.n_threads: number of threads for the current session
    % params.bb_tol: branch-and-bound tolerance (default 1e-4)
    % params.sdp_verbose: verbosity level of sdpnal+ (default 0)
    % params.sdp_tol: accuracy tolerance of sdpnal+ (default 1e-5)
    % params.cp_maxineq: maximum number of valid inequalities to separate (default 100000)
    % params.cp_maxiter: number of cutting-plane iterations (default 10)
    % params.cp_tol: cutting-plane tolerance (default 1e-4)
    % params.cp_percineq: percentage of inequalities to add (0.1 = 10%)
    % params.cp_epsineq: tolerance for checking the violation of inequalities (default 1e-4)
    % params.cp_activeineq: tolerance for checking active inequalities (default 1e-5) 
    % params.cp_inheritineq: inherit inequalitis from paret node (default yes - 1)
    % params.gurobi_vebose: verbosity level of gurobi (default 0)
    
    warning off;
    maxNumCompThreads(params.n_threads);
 
    result = struct();
    result.ub_list = [];
    result.lb_list = [];
    result.gap_list = [];
    result.time_list = [];

    [original_n, original_m] = size(W);
   
    [blk, At, b, C, L, shr_B_cell, T_U, ML_graph_U, T_V, ML_graph_V] = build_biclustering_shrinking(W, k, init_ML_U, init_CL_U, init_ML_V, init_CL_V, init_B_cell);
    n = size(T_U, 1); % new problem size
    m = size(T_V, 1); % new problem size
    result.n = n;
    result.m = m;
    
    Bt = svec(blk, shr_B_cell, 1);
    l = zeros(size(shr_B_cell, 2), 1);

    options.printlevel = params.sdp_verbose;
    options.tol = params.sdp_tol;
    options.stopoption = 0;
    options.AATsolve.method = 'direct';
    options.maxiter = 100000;
    options.maxtime = 7200;
    
    ub_max_eig = 2;
    [~, Yopt, ~, y, ~, Z2, y2, ~, info, ~] = sdpnalplus(blk, At, C, b, L, [], Bt, l, [], options);
    [~, ub] = safe_bound_error(blk, At, C, b, y, Z2, Bt, y2, l, ub_max_eig);
    fprintf('\n\n\t Upper bound = %10.9e \n', ub);
    
    Zshr = Yopt{1};
    Zuu = Zshr(1:n, 1:n);
    Zvv = Zshr(n+1:n+m, n+1:n+m); 
    Zuv = Zshr(1:n, n+1:n+m);
    Z = [T_U'*Zuu*T_U, T_U'*Zuv*T_V; T_V'*Zuv'*T_U, T_V'*Zvv*T_V];
    % heuristic with the SDP solution (provides a valid lower bound)
    [lb, Xu, Xv] = biclustering_heuristic(Z, W, k, params.gurobi_verbose);
    fprintf('\t Lower bound = %10.9e \n', lb);
         
    result.ub_list = [result.ub_list; ub];
    result.lb_list = [result.lb_list; lb];
    result.time_list = [result.time_list; info.totaltime];
    
    result.best_Z = Z;
    result.best_ub = ub;
    
    if lb < global_lb
        result.best_lb = global_lb;
        result.best_Xu = sparse(global_Xu);
        result.best_Xv = sparse(global_Xv);
    else
        result.best_lb = lb;
        result.best_Xu = sparse(Xu);
        result.best_Xv = sparse(Xv);
    end
    
    result.best_B_cell = init_B_cell;
    result.cp_iter = 0;
    result.cp_flag = -2;
    
    % -2 - maximum number of iterations
    % -1 - SDP not solved or partially solved successfully
    %  0 - no violated inequalities
    %  1 - node must be pruned
    %  2 - upper bound greater than the previous one
    %  3 - upper bound decrease is not sufficiently large
    
    gap = (result.best_ub - result.best_lb) / abs(result.best_ub);
    fprintf('\t Relative gap = %10.9e \n\n', gap);
    result.gap_list = [result.gap_list; gap];

    if gap <= params.bb_tol
        result.cp_flag = 1;
        return
    end

    % these inequalities should be inherited
    B_cell_or = init_B_cell;

    B_cell = shr_B_cell;
    l = zeros(size(shr_B_cell, 2), 1);
    
    for i=1:params.cp_maxiter
        
        % remove inactive inequalities
        if ~isempty(B_cell) && ~isempty(l)
            fprintf('\t Inequalities before removal = %d \n', size(B_cell, 2));
            Zvec = svec(blk, Zshr, 1);
            Bt = svec(blk, B_cell, 1);
            active_idx = abs(Bt{1}' * Zvec) <= params.cp_activeineq;
            B_cell = B_cell(active_idx);
            fprintf('\t Inequalities after removal = %d \n', size(B_cell, 2));
        end
        
        % separate pair and triangle inequalities
        [B_ineq_or, ~] = separate_inequalities(Z, original_n, original_m, params.cp_epsineq, params.cp_maxineq, params.cp_percineq);
        B_cell_or = [B_cell_or, B_ineq_or];

        B_ineq = shrink_cuts(ML_graph_U, ML_graph_V, B_ineq_or, n, m, original_n);
        n_ineq = size(B_ineq, 2);
        fprintf('\t Inequalities added = %d \n\n', n_ineq);
        
        if n_ineq <= n+m
            result.cp_flag = 0;
            break
        end
        
        B_cell = [B_cell, B_ineq];
        current_ineq = size(B_cell, 2);
        Bt = svec(blk, B_cell, 1);
        l = zeros(current_ineq, 1);
        % solve SDP with inequalities
        [~, Yopt, ~, y, ~, Z2, y2, ~, info, ~] = sdpnalplus(blk, At, C, b, L, [], Bt, l, [], options);
        result.cp_iter = result.cp_iter + 1;
        
        % SDP is not solved successfully
        if info.termcode == 1
            result.cp_flag = -1;
            break
        end
   
        [~, ub] = safe_bound_error(blk, At, C, b, y, Z2, Bt, y2, l, ub_max_eig);
        fprintf('\n\n\t Upper bound = %10.9e\n', ub);
        
        Zshr = Yopt{1};
        Zuu = Zshr(1:n, 1:n);
        Zvv = Zshr(n+1:n+m, n+1:n+m); 
        Zuv = Zshr(1:n, n+1:n+m);
        Z = [T_U'*Zuu*T_U, T_U'*Zuv*T_V; T_V'*Zuv'*T_U, T_V'*Zvv*T_V];
        % heuristic with the SDP solution (provides a lower bound)
        [lb, Xu, Xv] = biclustering_heuristic(Z, W, k, params.gurobi_verbose);
        fprintf('\t Lower bound = %10.9e\n', lb);
    
        result.ub_list = [result.ub_list; ub];
        result.lb_list = [result.lb_list; lb];
        result.time_list = [result.time_list; info.totaltime];
        
        if result.best_lb < lb
            % update best upper bound
            result.best_lb = lb;
            result.best_Xu = sparse(Xu); % best biclustering
            result.best_Xv = sparse(Xv); % best biclustering
        end
        
        % current upper bound greater than previous upper bound
        if ub > result.best_ub
            result.cp_flag = 2;
            gap = (result.best_ub - result.best_lb) / abs(result.best_ub);
            result.gap_list = [result.gap_list; gap];
            break
        end
        
        cp_stop = (result.best_ub - ub) / result.best_ub;
        % current bound is not sufficiently larger than previous lower bound
        if cp_stop <= params.cp_tol
            result.cp_flag = 3;
            result.best_ub = ub;
            result.best_Z = Z;
            result.best_B_cell = B_cell_or;
            gap = (result.best_ub - result.best_lb) / abs(result.best_ub);
            fprintf('\t Relative gap = %10.9e\n\n', gap);
            result.gap_list = [result.gap_list; gap];
            if gap <= params.bb_tol
                result.cp_flag = 1;
            end
            break
        end
        
        % update best upper bound
        result.best_ub = ub;
        result.best_Z = Z;
        result.best_B_cell = B_cell_or;
        
        gap = (result.best_ub - result.best_lb) / abs(result.best_ub);
        fprintf('\t Relative gap = %10.9e\n\n', gap);
        result.gap_list = [result.gap_list; gap];
        % prune the node
        if gap <= params.bb_tol
            result.cp_flag = 1;
            break
        end
                
    end
    

end
