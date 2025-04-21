function str = polyToHTMLString(coeffs, varSymbol)
    % POLYHTMLSTRING
    % --------------
    % This enhanced helper function converts a vector of polynomial coefficients 
    % into an HTML-formatted string. Exponents are rendered using <sup> tags.
    %
    % Input:
    %   coeffs - A numeric vector of polynomial coefficients.
    %   varSymbol - Optional. Variable symbol to use ('s' by default, could be 'z')
    % Output:
    %   str - An HTML-formatted string representing the polynomial.
    
    % Set default variable symbol if not provided
    if nargin < 2 || isempty(varSymbol)
        varSymbol = 's';
    end
    
    % Return 0 for empty or all-zero coefficients
    if isempty(coeffs) || all(coeffs == 0)
        str = '0';
        return;
    end
    
    % Polynomial order (highest power)
    n = length(coeffs) - 1;
    
    % Initialize array to store terms
    terms = {};
    
    % Process each coefficient
    for i = 1:length(coeffs)
        coeff = coeffs(i);
        exp = n - (i - 1);
        
        % Skip zero coefficients
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
        
        % Format coefficient (omit 1 for non-constant terms)
        if coeff == 1 && exp ~= 0
            coeffStr = '';
        else
            % Handle floating point numbers with precision
            if floor(coeff) == coeff
                coeffStr = num2str(coeff);
            else
                coeffStr = num2str(coeff, '%.4g');
            end
        end
        
        % Format the term with variable and exponent
        if exp > 1
            term = [signStr, coeffStr, varSymbol, '<sup>', num2str(exp), '</sup>'];
        elseif exp == 1
            term = [signStr, coeffStr, varSymbol];
        else
            term = [signStr, coeffStr];
        end
        
        terms{end+1} = term;
    end
    
    % Join all terms into a single string
    str = strjoin(terms, '');
    
    % If the first character is a space (from ' + '), remove it
    if ~isempty(str) && str(1) == ' '
        str = str(2:end);
    end
end