function [A_cell_Z_cl, b_cl] = add_cannot_link_Zuu(n, m, CL)
    
    n_cl = size(CL, 1);
    A_cell_Z_cl = cell(1, n_cl);
    b_cl = zeros(n_cl, 1);
    for c=1:n_cl
        i = CL(c, 1);
        j = CL(c, 2);
        A_cell_Z_cl{c} = sparse([i, j], [j, i], [0.5, 0.5], n+m, n+m);
    end
    
end