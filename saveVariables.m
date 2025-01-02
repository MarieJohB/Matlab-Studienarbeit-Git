function saveVariables(checkboxes)
    % Retrieve the values of the checkboxes
    selectedVars = {};
    for i = 1:numel(checkboxes)
        if checkboxes{i}.Value
            selectedVars{end+1} = checkboxes{i}.Text; %#ok<AGROW>
        end
    end