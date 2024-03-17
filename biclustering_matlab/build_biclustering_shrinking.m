function [blk, At, b, C, L, init_B_cell, T_U, ML_graph_U, T_V, ML_graph_V] = build_biclustering_shrinking(W, k, init_ML_U, init_CL_U, init_ML_V, init_CL_V, parent_B_cell)

    % W: weight matrix n x m
    % k: number of clusters
    % init_ML_U: matrix of must-link constraints between vertices in U (n_ml_U x 2)
    % init_CL_U: matrix of cannot-link constraints between vertices in U (n_cl_U x 2)
    % init_ML_V: matrix of must-link constraints between vertices in V (n_ml_V x 2)
    % init_CL_V: matrix of cannot-link constraints between vertices in V (n_cl_V x 2)
    % parent_Bcell: valid inequalities for the parent node
    % return problem data in SDPNAL+ format

    [original_n, original_m] = size(W);
    
    % transformation for shrinking U
    [T_U, ML_graph_U] = build_T(original_n, init_ML_U);
    % plot(ML_graph_U)
    n = size(T_U, 1); % new size
    CL_U = update_CL(ML_graph_U, init_CL_U);
    e_U = T_U * ones(original_n, 1);
    
    % transformation for shrinking V
    [T_V, ML_graph_V] = build_T(original_m, init_ML_V);
    % plot(ML_graph_V)
    m = size(T_V, 1); % new size
    CL_V = update_CL(ML_graph_V, init_CL_V);
    e_V = T_V * ones(original_m, 1);
    
    W_shr = T_U * W * T_V';
    W_full = 0.5*[zeros(n, n), W_shr; W_shr', zeros(m, m)];
    
    [Zuu_rowsum, Zuu_trace] = Z_slice_Zuu_shrinking(n, m, e_U);
    [Zvv_rowsum, Zvv_trace] = Z_slice_Zvv_shrinking(n, m, e_V);
    
    [A_cell_Z_cl_U, b_Z_cl_U] = add_cannot_link_Zuu(n, m, CL_U);
    [A_cell_Z_cl_V, b_Z_cl_V] = add_cannot_link_Zvv(n, m, CL_V);


    Acell = [Zuu_rowsum, Zuu_trace, Zvv_rowsum, Zvv_trace, A_cell_Z_cl_U, A_cell_Z_cl_V];
    b = [ones(n, 1); k; ones(m, 1); k; b_Z_cl_U; b_Z_cl_V];
        
    blk = cell(1);
    blk{1,1} = 's';
    blk{1,2} = n+m;
    
    At = svec(blk, Acell, 1);

    C = cell(1);
    C{1} = -sparse(W_full);

    L = cell(1);
    L{1} = 0; % we want Z >= 0
    
    % inherit cuts and change indices (note that in parent_B_cell the
    % indices are from 1 to original_n and from 1 to original_m, respectively)
    init_B_cell = shrink_cuts(ML_graph_U, ML_graph_V, parent_B_cell, n, m, original_n);
    
end