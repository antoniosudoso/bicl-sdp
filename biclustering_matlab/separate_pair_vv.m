function [B_pair, violations] = separate_pair_vv(Z, n, m, eps, max_pair)

    rng(1993);
    % m_idx = randperm(m);
    
    B_pair = cell(1, max_pair);
    violations = zeros(max_pair, 1);
    
    c = 1;
    stop = false;
    
    % PAIR INEQUALITIES Zvv
    
    for i=1:m
        for j=i+1:m
            viol1 = 0.5*Z(i+n,j+n) + 0.5*Z(j+n,i+n) - Z(i+n,i+n);
            viol2 = 0.5*Z(i+n,j+n) + 0.5*Z(j+n,i+n) - Z(j+n,j+n);
            if viol1 >= eps
                violations(c) = viol1;
                B_pair{c} = sparse([i+n; j+n; i+n], [j+n; i+n; i+n], [-0.5; -0.5; 1], n+m, n+m);
                c = c + 1;
            end
            if viol2 >= eps
                violations(c) = viol2;
                B_pair{c} = sparse([i+n; j+n; j+n], [j+n; i+n; j+n], [-0.5; -0.5; 1], n+m, n+m);
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