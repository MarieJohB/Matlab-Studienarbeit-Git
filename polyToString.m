function str = poly2str(p, var)
    % Check if p is zero or empty
    if isempty(p) || (length(p) == 1 && p(1) == 0)
        str = '0';
        return;
    end
    
    % Initialize an empty string
    str = '';
    
    % Loop through each coefficient
    n = length(p);
    for i = 1:n
        % Calculate the power for this coefficient
        power = n - i;
        
        % Skip zero coefficients (except for constant term if all zeros)
        if p(i) == 0 && (power > 0 || i < n)
            continue;
        end
        
        % Add the sign
        if ~isempty(str)
            if p(i) >= 0
                str = [str, ' + '];
            else
                str = [str, ' - '];
                p(i) = -p(i); % Make it positive for display
            end
        elseif p(i) < 0
            str = '-';
            p(i) = -p(i);
        end
        
        % Format the coefficient
        if power == 0
            % Constant term
            str = [str, num2str(p(i), '%.4g')];
        else
            % Term with variable
            if p(i) == 1
                % Coefficient of 1 is not shown (except for constant term)
                coef_str = '';
            else
                coef_str = num2str(p(i), '%.4g');
            end
            
            % Add the variable and exponent
            if power == 1
                str = [str, coef_str, var];
            else
                str = [str, coef_str, var, '^{', num2str(power), '}'];
            end
        end
    end
    
    % If the string is still empty (all zeros), return zero
    if isempty(str)
        str = '0';
    end
end