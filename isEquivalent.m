function result = isEquivalent(expr1, expr2)
    % Check if two expressions are equivalent
    % More robust than direct equality testing
    try
        % Try direct equality first
        result = (expr1 == expr2);
        
        % If that didn't work, try simplifying and comparing
        if ~result
            result = isequal(simplify(expr1 - expr2), 0);
        end
    catch
        % If all else fails, assume they're not equivalent
        result = false;
    end
end