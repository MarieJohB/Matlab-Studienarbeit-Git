% Helper function to format numbers consistently
function str = formatNumber(num, precision)
    if nargin < 2
        precision = 2; % Default precision
    end
    
    if isnan(num)
        str = 'N/A';
    elseif abs(num) > 1000 || (abs(num) < 0.01 && abs(num) > 0)
        % Use scientific notation for very large/small numbers
        str = sprintf(['%.' num2str(precision) 'e'], num);
    else
        % Use fixed-point notation with appropriate precision
        str = sprintf(['%.' num2str(precision) 'f'], num);
    end
end