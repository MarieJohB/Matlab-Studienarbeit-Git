function p = removeLeadingZeros(p)
    % Check if p is a cell and convert if needed
    if iscell(p)
        try
            p = cell2mat(p);
        catch
            % If conversion fails, just return as is
            return;
        end
    end
    
    % Find first non-zero element
    idx = 1;
    while idx <= length(p) && p(idx) == 0
        idx = idx + 1;
    end
    
    % Return appropriate result
    if idx > length(p)
        p = 0;  % All zeros
    else
        p = p(idx:end);  % Truncate leading zeros
    end
end