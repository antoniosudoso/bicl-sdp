function [B_cell, n_ineq] = separate_inequalities(Z, n, m, eps, max_sep, perc_ineq)

    B_cell = [];
    viol = [];

    [B_pair_uu, viol_pair_uu] = separate_pair_uu(Z, n, m, eps, max_sep);
    fprintf('\t Pair UU = %d \n', size(B_pair_uu, 2));
    B_cell = [B_cell, B_pair_uu];
    viol = [viol; viol_pair_uu];
    
    [B_pair_vv, viol_pair_vv] = separate_pair_vv(Z, n, m, eps, max_sep);
    fprintf('\t Pair VV = %d \n', size(B_pair_vv, 2));
    B_cell = [B_cell, B_pair_vv];
    viol = [viol; viol_pair_vv];
   
    [B_tri_uu, viol_tri_uu] = separate_triangle_uu(Z, n, m, eps, max_sep);
    fprintf('\t Triangle UU = %d \n', size(B_tri_uu, 2));
    B_cell = [B_cell, B_tri_uu];
    viol = [viol; viol_tri_uu];
    
    [B_tri_vv, viol_tri_vv] = separate_triangle_vv(Z, n, m, eps, max_sep);
    fprintf('\t Triangle VV = %d \n', size(B_tri_vv, 2));
    B_cell = [B_cell, B_tri_vv];
    viol = [viol; viol_tri_vv];
    
    %[B_pair_uv, viol_pair_uv] = separate_pair_uv(Z, n, m, eps, max_sep);
    %disp(size(B_pair_uv))
    %viol_pair_uv
    %B_cell = [B_cell, B_pair_uv];
    %viol = [viol; viol_pair_uv];
    
    %[B_tri_uv, viol_tri_uv] = separate_triangle_uv(Z, n, m, eps, max_sep);
    %disp(size(B_tri_uv))
    %viol_tri_uv
    %B_cell = [B_cell, B_tri_uv];
    %viol = [viol; viol_tri_uv];
    
    n_ineq = length(viol);
    max_ineq = floor(max_sep * perc_ineq);
        
    if n_ineq > max_ineq
        n_ineq = max_ineq;
        [~, id_sorted] = sort(viol, 'descend');
        %figure(1);
        %plot(viol_sorted);
        B_cell = B_cell(id_sorted);
        B_cell = B_cell(1:max_ineq);
    end
    
        
end