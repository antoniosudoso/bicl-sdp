function [B_pair, violations] = separate_pair_uu(Z, n, m, eps, max_pair)

    rng(1993);
    % n_idx = randperm(n);
    
    B_pair = cell(1, max_pair);
    violations = zeros(max_pair, 1);
    
    c = 1;
    stop = false;
    
    % PAIR INEQUALITIES Zuu
    
    for i=1:n
        for j=i+1:n
            viol1 = 0.5*Z(i,j) + 0.5*Z(j,i) - Z(i,i);
            viol2 = 0.5*Z(i,j) + 0.5*Z(j,i) - Z(j,j);
            if viol1 >= eps
                violations(c) = viol1;
                B_pair{c} = sparse([i; j; i], [j; i; i], [-0.5; -0.5; 1], n+m, n+m);
                c = c + 1;
            end
            if viol2 >= eps
                violations(c) = viol2;
                B_pair{c} = sparse([i; j; j], [j; i; j], [-0.5; -0.5; 1], n+m, n+m);
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