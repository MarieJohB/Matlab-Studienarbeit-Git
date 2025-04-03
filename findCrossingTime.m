% Utility function for finding rise/settling times
function t_cross = findCrossingTime(t, y, level)
    idx = find(y >= level, 1);
    if isempty(idx) || idx == 1
        t_cross = NaN;
    else
        % Linear interpolation for more accurate crossing time
        t_cross = interp1([y(idx-1), y(idx)], [t(idx-1), t(idx)], level);
    end
end