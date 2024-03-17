function [best_lb, best_Xu, best_Xv] = biclustering_heuristic(Z, W, k, verbose)

    [n, m] = size(W);
    
    %[l, s, r] = svd(Z);
    %Z = l(:, 1:k)*s(1:k, 1:k)*r(:, 1:k)';
                
    Zuu = Z(1:n, 1:n);
    Zvv = Z(n+1:n+m, n+1:n+m); 
    
    [labels_u, ~] = kmeans(Zuu, k, 'Replicates', 500, 'Start', 'plus');
    Xu = zeros(n, k);
    for i=1:n
        Xu(i, labels_u(i)) = 1;
    end
    
    Xu = Xu * diag(1./sqrt(sum(Xu, 1)));
    [labels_v, ~] = kmeans(Zvv, k, 'Replicates', 500, 'Start', 'plus');
    Xv = zeros(m, k);
    for i=1:m
        Xv(i, labels_v(i)) = 1;
    end
    Xv = Xv * diag(1./sqrt(sum(Xv, 1)));
    
    % Index helper function
    assidx = @(i, j) i+(j-1)*k;
    
    % Build model
    model.A = sparse(2*k, k*k);
    model.lb = zeros(k*k, 1);
    obj = zeros(k, k);

    for i=1:k
        for j=1:k
            model.A(i, assidx(i, j)) = 1;
            obj(i, j) = Xu(:, i)'*W*Xv(:, j);
        end
    end
        
    %disp(obj)
    
    %obj = obj';
        
    for j=1:k
        model.A(j+k, 1+k*(j-1):k*j) = 1;
    end
        
    model.rhs = ones(2*k, 1);
    model.sense = repmat('=', 2*k, 1);
            
    model.modelsense = 'max';
    model.obj = obj(:);
    params.outputflag = verbose;
    result = gurobi(model, params);
    %disp(result.objval);
    delta = reshape(result.x, [k, k]);
    %disp(delta)
    [r_idx, c_idx] = find(delta);
        
    % recover best normalized partition matrices
    best_Xu = zeros(n, k);
    best_Xv = zeros(m, k);
    for i=1:k
        best_Xu(:, i) = Xu(:, r_idx(i));
        best_Xv(:, i) = Xv(:, c_idx(i));   
    end
    
    % old_lb = trace(Xu'*W*Xv);
    best_lb = trace(best_Xu'*W*best_Xv);
     
end