function str = formatNumberAsFraction(num)
    % Format a number as a fraction if it's a decimal, or as a regular number otherwise
    if isreal(num) && mod(num, 1) ~= 0
        % Convert decimal to fraction
        [n, d] = rat(num, 1e-6);
        if d > 1 && d <= 100  % Reasonable denominator limit
            str = [num2str(n) '/' num2str(d)];
        else
            str = num2str(num);
        end
    else
        str = num2str(num);
    end
end