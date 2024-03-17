function [branching_type, i_idx, j_idx, max_val] = get_branching_pair(Z, n, m)

    Zuu = Z(1:n, 1:n);
    Zvv = Z(n+1:n+m, n+1:n+m);
    
    % [i, j, min_val] uu block
    values_U = zeros(n*(n-1)/2, 3);
    count = 1;
    for i=1:n
        for j=i+1:n
            min_val = min(Zuu(i, j), Zuu(i, i) - Zuu(i, j));
            values_U(count, :) = [i, j, min_val];
            count = count + 1;
        end
    end
     
    % [i, j, min_val] vv block
    values_V = zeros(m*(m-1)/2, 3);
    count = 1;
    for i=1:m
        for j=i+1:m
            min_val = min(Zvv(i, j), Zvv(i, i) - Zvv(i, j));
            values_V(count, :) = [i, j, min_val];
            count = count + 1;
        end
    end
    
    sorted_values_Zuu = sortrows(values_U, 3, 'descend');
    sorted_values_Zvv = sortrows(values_V, 3, 'descend');
    max_val_U = sorted_values_Zuu(1, 3);
    max_val_V = sorted_values_Zvv(1, 3);
        
    if n*max_val_U > m*max_val_V 
        % branch on U
        branching_type = 0;
        max_val = max_val_U;
        i_idx = sorted_values_Zuu(1, 1);
        j_idx = sorted_values_Zuu(1, 2);
    else
        % branch on V
        branching_type = 1;
        max_val = max_val_V;
        i_idx = sorted_values_Zvv(1, 1);
        j_idx = sorted_values_Zvv(1, 2);
    end
    
    if j_idx < i_idx
        temp_idx = j_idx;
        j_idx = i_idx;
        i_idx = temp_idx;
    end
    
    if max_val_U < 1e-4 && max_val_V < 1e-4
        i_idx = -1;
        j_idx = -1;
    end

end