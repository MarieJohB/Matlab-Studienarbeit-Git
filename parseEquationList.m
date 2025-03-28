function eqList = parseEquationList(eqListStr)
    % Parse a semicolon-separated list of equations with improved robustness
    
    if isempty(eqListStr)
        eqList = {};
        return;
    end
    
    % Split by semicolons
    items = split(eqListStr, ';');
    
    % Trim whitespace and filter empty items
    eqList = {};
    for i = 1:length(items)
        item = strtrim(items{i});
        if ~isempty(item)
            eqList{end+1} = item;
        end
    end
end