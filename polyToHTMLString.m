function str = polyToHTMLString(coeffs)
    % polyToHTMLString
    % -----------------
    % This helper function converts a vector of polynomial coefficients into an
    % HTML-formatted string. Exponents are rendered using <sup> tags.
    %
    % Input:
    %   coeffs - A numeric vector of polynomial coefficients.
    % Output:
    %   str    - An HTML-formatted string representing the polynomial.
    
    if isempty(coeffs) || all(coeffs == 0)
        str = '0';
        return;
    end
    
    n = length(coeffs) - 1;
    terms = {};
    for i = 1:length(coeffs)
        coeff = coeffs(i);
        exp = n - (i - 1);
        if coeff == 0
            continue;
        end
        % Determine sign
        if coeff < 0
            signStr = ' - ';
            coeff = abs(coeff);
        else
            if isempty(terms)
                signStr = '';
            else
                signStr = ' + ';
            end
        end
        % Omit coefficient 1 for non-constant terms
        if coeff == 1 && exp ~= 0
            coeffStr = '';
        else
            coeffStr = num2str(coeff);
        end
        if exp > 1
            term = [signStr, coeffStr, 's<sup>', num2str(exp), '</sup>'];
        elseif exp == 1
            term = [signStr, coeffStr, 's'];
        else
            term = [signStr, coeffStr];
        end
        terms{end+1} = term;
    end
    str = strjoin(terms, '');
end