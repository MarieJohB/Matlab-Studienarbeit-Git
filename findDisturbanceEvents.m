% Helper function to find disturbance event times
function event_times = findDisturbanceEvents(t, d_signal)
    % Find significant changes in the disturbance signal
    d_diff = abs(diff(d_signal));
    threshold = max(d_diff) * 0.1; % 10% of max change is significant
    
    % Find indices where significant changes occur
    event_indices = find(d_diff > threshold);
    
    % Convert indices to times
    event_times = t(event_indices);
    
    % If too many events, limit to the most significant ones
    if length(event_times) > 5
        [~, sorted_idx] = sort(d_diff(event_indices), 'descend');
        event_times = event_times(sorted_idx(1:5));
        event_times = sort(event_times); % Re-sort by time
    end
end