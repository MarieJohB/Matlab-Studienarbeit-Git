function arg = extractFunctionArg(expr, funcName)
    % Extract the argument from a function call like sin, cos, exp, etc.
    % This helps get the actual parameter that might be a complex expression
    
    % Convert to string and use regex to extract the argument
    exprStr = char(expr);
    
    % For common transcendental functions
    pattern = [funcName, '\((.*)\)'];
    matches = regexp(exprStr, pattern, 'tokens');
    
    if ~isempty(matches) && ~isempty(matches{1})
        % Convert the extracted argument back to symbolic
        arg = sym(matches{1}{1});
    else
        % Fallback to just using symvar if regex extraction fails
        arg = symvar(expr);
    end
end