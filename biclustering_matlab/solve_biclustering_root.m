function [result] = solve_biclustering_root(W, k, params)

    % SOLVE SDP RELAXATION FOR BICLUSTERING
    
    % W: weight matrix (n x m)
    % k: number of clusters
    
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

    %disp(params);
    
    warning off;    
    maxNumCompThreads(params.n_threads);
 
    result = struct();
    result.ub_list = [];
    result.lb_list = [];
    result.gap_list = [];
    result.time_list = [];

    [n, m] = size(W);
    [blk, At, b, C, L] = build_biclustering(W, k);
    result.n = n;
    result.m = m;
        
    options.printlevel = params.sdp_verbose;
    options.tol = params.sdp_tol;
    options.stopoption = 0;
    options.AATsolve.method = 'direct';
    options.maxiter = 100000;
    options.maxtime = 7200;
    
    ub_max_eig = 2;
    [~, Yopt, ~, y, ~, Z2, ~, ~, info, ~] = sdpnalplus(blk, At, C, b, L, [], [], [], [], options);
    [~, ub] = safe_bound_error(blk, At, C, b, y, Z2, [], 0, 0, ub_max_eig);
    fprintf('\n\n\t Upper bound = %10.9e \n', ub);
    
    Z = Yopt{1};
    % heuristic with the SDP solution (provides a lower bound)
    [lb, Xu, Xv] = biclustering_heuristic(Z, W, k, params.gurobi_verbose);
    fprintf('\t Lower bound = %10.9e \n', lb);
         
    result.ub_list = [result.ub_list; ub];
    result.lb_list = [result.lb_list; lb];
    result.time_list = [result.time_list; info.totaltime];
    
    result.best_Z = Z;
    result.best_ub = ub;
    result.best_lb = lb;
    result.best_Xu = sparse(Xu); % best biclustering
    result.best_Xv = sparse(Xv); % best biclustering
    result.best_B_cell = cell(0);
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

    B_cell = [];
    l = [];
    
    for i=1:params.cp_maxiter
        
        % remove inactive inequalities
        if ~isempty(B_cell) && ~isempty(l)
            fprintf('\t Inequalities before removal = %d \n', size(B_cell, 2));
            Zvec = svec(blk, Z, 1);
            Bt = svec(blk, B_cell, 1);
            active_idx = abs(Bt{1}' * Zvec) <= params.cp_activeineq;
            B_cell = B_cell(active_idx);
            fprintf('\t Inequalities after removal = %d \n', size(B_cell, 2));
        end
        
        [B_ineq, n_ineq] = separate_inequalities(Z, n, m, params.cp_epsineq, params.cp_maxineq, params.cp_percineq);
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
        
        Z = Yopt{1};
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
            result.best_B_cell = B_cell;
            gap = (result.best_ub - result.best_lb) / abs(result.best_ub);
            fprintf('\t Relative gap = %10.9e\n\n', gap);
            result.gap_list = [result.gap_list; gap];
            if gap <= params.bb_tol
                result.cp_flag = 1;
                break
            end
            break
        end
        
        % update best upper bound
        result.best_ub = ub;
        result.best_Z = Z;
        result.best_B_cell = B_cell;
        
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
