function varList = parseVariableList(varListStr)
    % Parse a comma-separated list of variables with improved robustness
    
    if isempty(varListStr)
        varList = {};
        return;
    end
    
    % Split by commas
    items = split(varListStr, ',');
    
    % Trim whitespace and filter empty items
    varList = {};
    for i = 1:length(items)
        item = strtrim(items{i});
        if ~isempty(item)
            % Check if valid variable name
            if isvarname(item)
                varList{end+1} = item;
            else
                warning(['Invalid variable name skipped: ' item]);
            end
        end
    end
end