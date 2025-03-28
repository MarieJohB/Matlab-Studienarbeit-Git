function str = matrix2str(A)
    % Convert a matrix to a string representation
    % Format: [a b c; d e f] becomes "a b c; d e f"
    
    [m, n] = size(A);
    str = '';
    
    for i = 1:m
        row_str = '';
        for j = 1:n
            if j > 1
                row_str = [row_str, ' '];
            end
            row_str = [row_str, num2str(A(i,j))];
        end
        
        if i < m
            row_str = [row_str, '; '];
        end
        
        str = [str, row_str];
    end
end