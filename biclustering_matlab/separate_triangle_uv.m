function [B_triangle, violations] = separate_triangle_uv(Z, n, m, eps, max_triangle)
    
    rng(1993);
    n_idx = randperm(n);
    m_idx = randperm(m);
    
    B_triangle = cell(1, max_triangle);    
    violations = zeros(max_triangle, 1);
    
    c = 1;
    stop = false;
    
    % TRIANGLE INEQUALITIES Zuu, Zuv 
    
    for i=n_idx
        for j=n_idx
            if (i ~= j)
                for t=m_idx
                    viol1 = 0.5*Z(i,t+n) + 0.5*Z(t+n,i) + 0.5*Z(i,j) + 0.5*Z(j,i) - Z(i,i) - 0.5*Z(j,t+n) - 0.5*Z(t+n,j);
                    viol2 = 0.5*Z(i,t+n) + 0.5*Z(t+n,i) + 0.5*Z(i,j) + 0.5*Z(j,i) - Z(j,j) - 0.5*Z(j,t+n) - 0.5*Z(t+n,j);
                    if viol1 >= eps
                        violations(c) = viol1;
                        B_triangle{c} = sparse([i, t+n, i, j, i, j, t+n], [t+n, i, j, i, i, t+n, j], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+m, n+m);
                        c = c + 1;
                     end
                     if viol2 >= eps
                        violations(c) = viol2;
                        B_triangle{c} = sparse([i, t+n, i, j, j, j, t+n], [t+n, i, j, i, j, t+n, j], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+m, n+m);
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
    
    % TRIANGLE INEQUALITIES Zuu, Zuv 
    
    c = 1;
    stop = false;
    
    for i=m_idx
        for j=m_idx
            if (i ~= j)
                for t=n_idx
                    viol1 = 0.5*Z(i+n,t) + 0.5*Z(t,i+n) + 0.5*Z(i+n,j+n) + 0.5*Z(j+n,i+n) - Z(i+n,i+n) - 0.5*Z(j+n,t) - 0.5*Z(t,j+n);
                    viol2 = 0.5*Z(i+n,t) + 0.5*Z(t,i+n) + 0.5*Z(i+n,j+n) + 0.5*Z(j+n,i+n) - Z(j+n,j+n) - 0.5*Z(j+n,t) - 0.5*Z(t,j+n);
                    if viol1 >= eps
                        violations(c) = viol1;
                        B_triangle{c} = sparse([i+n, t, i+n, j+n, i+n, j+n, t], [t, i+n, j+n, i+n, i+n, t, j+n], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+m, n+m);
                        c = c + 1;
                     end
                     if viol2 >= eps
                        violations(c) = viol2;
                        B_triangle{c} = sparse([i+n, t, i+n, j+n, j+n, j+n, t], [t, i+n, j+n, i+n, j+n, t, j+n], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+m, n+m);
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