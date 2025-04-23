function fracsym = convertToFractions(sym_expr)
    % Convert decimal coefficients in a symbolic expression to fractions
    % This recursive function handles nested expressions
    
    % For atomic symbolic expressions (numbers, variables, etc.)
    if isempty(symvar(sym_expr)) && ~isa(sym_expr, 'symfun')
        % Check if it's a decimal number
        num_val = double(sym_expr);
        if isreal(num_val) && mod(num_val, 1) ~= 0
            % Convert decimal to fraction with tolerance
            [n, d] = rat(num_val, 1e-6);
            if d > 1 && d <= 100  % Only use fractions with reasonable denominators
                fracsym = sym(n)/sym(d);
            else
                fracsym = sym_expr;
            end
        else
            fracsym = sym_expr;
        end
        return;
    end
    
    % Special case for expressions with e^(a*t) format which are common in control systems
    try
        if has(sym_expr, exp(sym('t')*sym('a'))) || ...       % Check for e^(a*t)
           has(sym_expr, exp(sym('a')*sym('t')))              % Same but different order
            % Try direct string manipulation for better handling of exponentials
            exprStr = char(sym_expr);
            % Find all decimal numbers in the expression
            [starts, ends, tokens] = regexp(exprStr, '(\d+\.\d+)', 'start', 'end', 'tokens');
            
            % Process each decimal number
            offsetAdjustment = 0;
            for i = 1:length(starts)
                decimalValue = str2double(tokens{i}{1});
                [n, d] = rat(decimalValue, 1e-6);
                
                % Only replace if the denominator is reasonable
                if d > 1 && d <= 100
                    originalLength = ends(i) - starts(i) + 1;
                    fractionStr = [num2str(n) '/' num2str(d)];
                    
                    % Adjust for the new string length
                    adjustedStart = starts(i) + offsetAdjustment;
                    adjustedEnd = ends(i) + offsetAdjustment;
                    
                    % Replace the decimal with the fraction in the string
                    exprStr = [exprStr(1:adjustedStart-1), fractionStr, exprStr(adjustedEnd+1:end)];
                    
                    % Update the offset adjustment
                    offsetAdjustment = offsetAdjustment + (length(fractionStr) - originalLength);
                end
            end
            
            % Convert back to symbolic
            fracsym = sym(exprStr);
            return;
        end
    catch
        % If the special case handling fails, continue with the general approach
    end
    
    % Handle different types of symbolic expressions
    [op, args] = extractOperationAndArgs(sym_expr);
    
    if isempty(op)
        % It's a variable or another simple expression
        fracsym = sym_expr;
    else
        % Process each argument recursively
        for i = 1:length(args)
            args{i} = convertToFractions(args{i});
        end
        
        % Reconstruct the expression
        fracsym = reconstructExpression(op, args);
    end
end