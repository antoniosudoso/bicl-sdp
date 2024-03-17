function [B_triangle, violations] = separate_triangle_uu(Z, n, m, eps, max_triangle)
    
    rng(1993);
    n_idx = randperm(n);
    
    B_triangle = cell(1, max_triangle);    
    violations = zeros(max_triangle, 1);
    
    c = 1;
    stop = false;
    
    % TRIANGLE INEQUALITIES Zuu
    
    for i=n_idx
        for j=i+1:n
            for t=j+1:n
                 if (i ~= t)
                     viol1 = 0.5*Z(i,j) + 0.5*Z(j,i) + 0.5*Z(i,t) + 0.5*Z(t,i) - Z(i,i) - 0.5*Z(j,t) - 0.5*Z(t,j);
                     viol2 = 0.5*Z(i,j) + 0.5*Z(j,i) + 0.5*Z(j,t) + 0.5*Z(t,j) - Z(j,j) - 0.5*Z(i,t) - 0.5*Z(t,i);
                     viol3 = 0.5*Z(i,t) + 0.5*Z(t,i) + 0.5*Z(j,t) + 0.5*Z(t,j) - Z(t,t) - 0.5*Z(i,j) - 0.5*Z(j,i);
                     if viol1 >= eps
                        violations(c) = viol1;
                        B_triangle{c} = sparse([i, j, i, t, i, j, t], [j, i, t, i, i, t, j], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+m, n+m);
                        c = c + 1;
                     end
                     if viol2 >= eps
                        violations(c) = viol2;
                        B_triangle{c} = sparse([i, j, j, t, j, i, t], [j, i, t, j, j, t, i], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+m, n+m);
                        c = c + 1;
                     end
                     if viol3 >= eps
                        violations(c) = viol3;
                        B_triangle{c} = sparse([i, t, j, t, t, i, j], [t, i, t, j, t, j, i], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+m, n+m);
                        c = c + 1;
                     end
                     
                     if c-1 == max_triangle 
                         stop = true;
                         break;
                     end   
                 end
            end
            if stop == true
                break;
            end
        end
        if stop == true
            break;
        end
    end
    
    B_triangle = B_triangle(1:(c-1));
    violations = violations(1:(c-1));
    
end