function new_Bcell = shrink_cuts(ML_graph_U, ML_graph_V, init_B_cell, n, m, original_n)
    
    n_ineq = size(init_B_cell, 2);
    
    new_Bcell = cell(1, n_ineq);
    [bins_v_U , ~] = conncomp(ML_graph_U);
    [bins_v_V , ~] = conncomp(ML_graph_V);

    counter = 1;
    for c=1:n_ineq
        
        [id_i, id_j, v] = find(init_B_cell{c});
                
        d = size(id_i, 1);
        
        % these indices must be updated according to ML_graph_U
        if sum(id_i <= original_n) == d
            for t=1:d
                id_i(t) = bins_v_U(id_i(t));
                id_j(t) = bins_v_U(id_j(t));
            end
            new_Bcell{counter} = sparse(id_i, id_j, v, n+m, n+m);
            counter = counter + 1;
        end
        
        % these indices must be adapted according to ML_graph_V
        if sum(id_i >= original_n + 1) == d
            for t=1:d
                id_i(t) = bins_v_V(id_i(t)-original_n);
                id_j(t) = bins_v_V(id_j(t)-original_n);
            end
            new_Bcell{counter} = sparse(id_i+n, id_j+n, v, n+m, n+m);
            counter = counter + 1;
        end
        
                
    end
    
    n_ineq = counter - 1;
    new_Bcell = new_Bcell(1:n_ineq);
    
end