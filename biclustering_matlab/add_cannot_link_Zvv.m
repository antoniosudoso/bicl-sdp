function [A_cell_Z_cl, b_cl] = add_cannot_link_Zvv(n, m, CL)
    
    n_cl = size(CL, 1);
    A_cell_Z_cl = cell(1, n_cl);
    b_cl = zeros(n_cl, 1);
    for c=1:n_cl
        i = CL(c, 1);
        j = CL(c, 2);
        A_cell_Z_cl{c} = sparse([i + n, j + n], [j + n, i + n], [0.5, 0.5], n+m, n+m);
    end
    
end