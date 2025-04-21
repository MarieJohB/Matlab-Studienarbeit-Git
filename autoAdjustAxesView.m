function autoAdjustAxesView(ax, t, y)
    % autoAdjustAxesView
    % -----------------
    % Helper function to automatically adjust the axes view to focus on the response curve
    %
    % Parameters:
    % - ax: Axes handle to adjust
    % - t: Time vector of the response
    % - y: Amplitude vector of the response
    
    % Get the data range with some padding
    yMin = min(y);
    yMax = max(y);
    
    % Add padding (10% of the data range)
    yRange = yMax - yMin;
    if yRange < eps
        % If range is too small (constant function), use default padding
        padding = 0.1;
        yMin = yMin - padding;
        yMax = yMax + padding;
    else
        padding = 0.1 * yRange;
        yMin = yMin - padding;
        yMax = yMax + padding;
    end
    
    % Find the settling time (where the response stays within 2% of final value)
    finalValue = y(end);
    settlingBand = 0.02 * abs(finalValue);
    
    % Find the last time the response leaves the settling band
    settleTimes = find(abs(y - finalValue) > settlingBand);
    
    % Set the x-axis limit based on settling time or the entire response
    if ~isempty(settleTimes) && max(settleTimes) < length(t)
        tSettle = t(max(settleTimes) + 1);
        % Add 20% more time to show the settled response
        tMax = tSettle * 1.2;
        % Ensure we don't exceed the time vector
        tMax = min(tMax, t(end));
    else
        % If no settling time found or it's at the end, show the entire response
        tMax = t(end);
    end
    
    % Set the axes limits
    xlim(ax, [0, tMax]);
    ylim(ax, [yMin, yMax]);
    
    % Refresh the plot
    drawnow;
end