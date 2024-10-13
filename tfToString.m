function [num_str, den_str] = tfToString(num, den)
    % Convert the numerator to a string
    num_str = polyToString(num);
    
    % Convert the denominator to a string
    den_str = polyToString(den);
end