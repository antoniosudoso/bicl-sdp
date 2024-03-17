function [B_triangle, violations] = separate_triangle_vv(Z, n, m, eps, max_triangle)
    
    rng(1993);
    m_idx = randperm(m);
    
    B_triangle = cell(1, max_triangle);    
    violations = zeros(max_triangle, 1);
    
    c = 1;
    stop = false;
    
    % TRIANGLE INEQUALITIES Zvv
    
    for i=m_idx
        for j=i+1:m
            for t=j+1:m
                 if (i ~= t)
                     viol1 = 0.5*Z(i+n,j+n) + 0.5*Z(j+n,i+n) + 0.5*Z(i+n,t+n) + 0.5*Z(t+n,i+n) - Z(i+n,i+n) - 0.5*Z(j+n,t+n) - 0.5*Z(t+n,j+n);
                     viol2 = 0.5*Z(i+n,j+n) + 0.5*Z(j+n,i+n) + 0.5*Z(j+n,t+n) + 0.5*Z(t+n,j+n) - Z(j+n,j+n) - 0.5*Z(i+n,t+n) - 0.5*Z(t+n,i+n);
                     viol3 = 0.5*Z(i+n,t+n) + 0.5*Z(t+n,i+n) + 0.5*Z(j+n,t+n) + 0.5*Z(t+n,j+n) - Z(t+n,t+n) - 0.5*Z(i+n,j+n) - 0.5*Z(j+n,i+n);
                     if viol1 >= eps
                        violations(c) = viol1;
                        B_triangle{c} = sparse([i+n, j+n, i+n, t+n, i+n, j+n, t+n], [j+n, i+n, t+n, i+n, i+n, t+n, j+n], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+m, n+m);
                        c = c + 1;
                     end
                     if viol2 >= eps
                        violations(c) = viol2;
                        B_triangle{c} = sparse([i+n, j+n, j+n, t+n, j+n, i+n, t+n], [j+n, i+n, t+n, j+n, j+n, t+n, i+n], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+m, n+m);
                        c = c + 1;
                     end
                     if viol3 >= eps
                        violations(c) = viol3;
                        B_triangle{c} = sparse([i+n, t+n, j+n, t+n, t+n, i+n, j+n], [t+n, i+n, t+n, j+n, t+n, j+n, i+n], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+m, n+m);
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