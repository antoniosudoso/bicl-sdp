function [T, G] = build_T(n, ML)
    
    n_ml = size(ML, 1);
    G = graph();
    G = addnode(G, n);
    for c=1:n_ml
        G = addedge(G, ML(c, 1), ML(c, 2));
    end
    %plot(G);
    [bins, binsizes] = conncomp(G, 'OutputForm','cell');
    % bins is a cell array, and bins{j} contains the node IDs for all nodes that belong to component j.
    m = size(binsizes, 2);
    T = zeros(m, n);
    for i=1:m
        for j=bins{i}
            T(i, j) = 1;
        end
    end

end