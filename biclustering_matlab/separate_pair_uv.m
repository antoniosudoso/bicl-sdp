function [B_pair, violations] = separate_pair_uv(Z, n, m, eps, max_pair)

    rng(1993);
    % n_idx = randperm(n);
    % m_idx = randperm(m);
    
    B_pair = cell(1, max_pair);
    violations = zeros(max_pair, 1);
    
    c = 1;
    stop = false;
    
    % PAIR INEQUALITIES Zuu and Zuv
    
    for i=1:n
        for j=1:m
            viol = 0.5*Z(i,j+n) + 0.5*Z(j+n,i) - Z(i,i);
            if viol >= eps
                violations(c) = viol;
                B_pair{c} = sparse([i; j+n; i], [j+n; i; i], [-0.5; -0.5; 1], n+m, n+m);
                c = c + 1;
            end
            if c-1 == max_pair
                stop = true;
                break;
            end
        end
        if stop == true
            break;
        end
    end
    
    stop = false;
    
    for i=1:m
        for j=1:n
            viol = 0.5*Z(i+n,j) + 0.5*Z(j,i+n) - Z(i+n,i+n);
            if viol >= eps
                violations(c) = viol;
                B_pair{c} = sparse([i+n; j; i+n], [j; i+n; i+n], [-0.5; -0.5; 1], n+m, n+m);
                c = c + 1;
            end
            if c-1 == max_pair
                stop = true;
                break;
            end
        end
        if stop == true
            break;
        end
    end
    
    B_pair = B_pair(1:(c-1));
    violations = violations(1:(c-1));
    
        
end