function [A_rowsum, A_trace] = Z_slice_Zuu(n, m)

    A_rowsum = cell(1, n);
    A_trace = cell(1, 1);
    A_trace{1} = sparse([eye(n), zeros(n, m); zeros(m, n), zeros(m, m)]);
    
    %disp(full(A_trace{1}))
    
    c = 1;
    for i = 1:n
        A_rowsum{c} = sparse([repelem(i, n), 1:n], [1:n, repelem(i, n)], 0.5, n+m, n+m);
        %disp(full(A_rowsum{c}))
        c = c + 1;
    end
    
end