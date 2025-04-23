function rounded_expr = roundSymbolicNumerics(expr)
    % ROUNDSYMBOLICNUMERICS
    % ---------------------
    % Helper function to round all numeric values in a symbolic expression
    % to 5 decimal places, while preserving mathematical functions and variables.
    %
    % Input:
    %   expr - Symbolic expression
    % Output:
    %   rounded_expr - Expression with rounded numeric values
    
    % Process differently based on expression type
    if isa(expr, 'sym') && ~isAllVars(expr) && ~isempty(symvar(expr))
        % Get children of the expression
        [op, args] = children(expr);
        
        % Recursively process each argument
        for i = 1:length(args)
            args{i} = roundSymbolicNumerics(args{i});
        end
        
        % Reconstruct the expression
        rounded_expr = op(args{:});
    elseif isa(expr, 'sym') && isAllVars(expr)
        % It's a variable (like t), leave it as is
        rounded_expr = expr;
    elseif isa(expr, 'sym') && isempty(symvar(expr))
        % It's a pure number in symbolic form
        % Convert to double, round, and convert back to symbolic
        numeric_value = double(expr);
        if isinf(numeric_value) || isnan(numeric_value)
            % Special cases: infinity or NaN
            rounded_expr = expr;
        else
            % Round regular numbers
            rounded_value = round(numeric_value * 10^5) / 10^5;
            rounded_expr = sym(rounded_value);
        end
    else
        % For any other case, leave as is
        rounded_expr = expr;
    end
end