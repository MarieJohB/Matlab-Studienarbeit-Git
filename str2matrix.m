function A = str2matrix(str)
    % Convert a string to a matrix with enhanced error handling
    % Format: "a b c; d e f" becomes [a b c; d e f]
    
    % Handle empty input case
    if isempty(str)
        A = [];
        return;
    end
    
    % Handle different input types
    if iscell(str)
        % Handle cell array input (from text area with multiple lines)
        combined_str = '';
        for i = 1:length(str)
            if i > 1
                combined_str = [combined_str, ';'];
            end
            combined_str = [combined_str, str{i}];
        end
        str = combined_str;
    elseif isstring(str)
        str = char(str);
    elseif ~ischar(str)
        error('Input must be either a character vector, string scalar, or cell array of strings.');
    end
    
    % Replace commas with spaces if present
    str = strrep(str, ',', ' ');
    
    % Split by semicolons to get rows
    rows = strsplit(strtrim(str), ';');
    
    % Process each row
    A = [];
    for i = 1:length(rows)
        % Split by spaces to get columns
        row_str = strtrim(rows{i});
        if isempty(row_str)
            continue; % Skip empty rows
        end
        
        % Use regexp to handle multiple spaces
        elements = regexp(row_str, '\s+', 'split');
        
        % Filter out empty elements
        elements = elements(~cellfun(@isempty, elements));
        
        if isempty(elements)
            continue; % Skip if no valid elements
        end
        
        % Convert to double
        row_values = zeros(1, length(elements));
        for j = 1:length(elements)
            row_values(j) = str2double(elements{j});
            
            % Check for conversion errors
            if isnan(row_values(j)) && ~strcmp(elements{j}, 'NaN')
                warning('Invalid numeric value: %s in row %d, column %d. Using 0.', elements{j}, i, j);
                row_values(j) = 0;
            end
        end
        
        % Check if this row would make the matrix non-rectangular
        if ~isempty(A) && size(A, 2) ~= length(row_values)
            error('Matrix rows must have the same number of columns.');
        end
        
        % Add to matrix
        A = [A; row_values];
    end
    
    % Check if we have an empty matrix
    if isempty(A)
        % Return an explicitly empty matrix rather than []
        A = zeros(0, 0);
    end
end