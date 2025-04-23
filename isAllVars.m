function result = isAllVars(expr)
    % Check if expression is just a variable
    result = false;
    try
        result = (numel(symvar(expr)) == 1) && (char(expr) == char(symvar(expr)));
    catch
        % If there's an error, it's not a simple variable
        result = false;
    end
end