function str = coeffsToString(coeffs)
% Convert coefficient array to string
str = '';
for i = 1:length(coeffs)
    if i > 1
        str = [str ' + '];
    end
    
    if length(coeffs) - i > 0
        if length(coeffs) - i > 1
            str = [str sprintf('%.4g·s^%d', coeffs(i), length(coeffs)-i)];
        else
            str = [str sprintf('%.4g·s', coeffs(i))];
        end
    else
        str = [str sprintf('%.4g', coeffs(i))];
    end
end
end