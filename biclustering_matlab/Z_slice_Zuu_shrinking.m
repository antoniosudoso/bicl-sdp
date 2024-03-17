function [A_rowsum, A_trace] = Z_slice_Zuu_shrinking(n, m, e)

    % e_U = T_U 1_originaln (new size is n) Diag(e_U) = T_U T_U^\top

    A_rowsum = cell(1, n);
    A_trace = cell(1, 1);
    A_trace{1} = sparse([diag(e), zeros(n, m); zeros(m, n), zeros(m, m)]);
    
    %disp(full(A_trace{1}))
    
    c = 1;
    for i = 1:n
        A_rowsum{c} = sparse([repelem(i, n), 1:n], [1:n, repelem(i, n)], [0.5 .* ones(n, 1) .* e; 0.5 .* ones(n, 1) .* e], n+m, n+m);
        %disp(full(A_rowsum{c}))
        c = c + 1;
    end
    
end