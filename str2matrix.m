function A = str2matrix(str)
    % Convert a string to a matrix
    % Format: "a b c; d e f" becomes [a b c; d e f]
    
    % Replace commas with spaces if present
    str = strrep(str, ',', ' ');
    
    % Split by semicolons to get rows
    rows = strtrim(split(str, ';'));
    
    % Process each row
    A = [];
    for i = 1:length(rows)
        % Split by spaces to get columns
        elements = strtrim(split(rows{i}));
        
        % Convert to double
        row_values = zeros(1, length(elements));
        for j = 1:length(elements)
            row_values(j) = str2double(elements{j});
        end
        
        % Add to matrix
        A = [A; row_values];
    end
end