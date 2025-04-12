function createStabilityTab(tab, batchResults)
% Create axes for stability plot
stabAxes = uiaxes(tab, 'Position', [50 350 900 300]);

% Plot stability
paramValues = batchResults.paramValues;
stabilityValues = double(batchResults.stability);

plot(stabAxes, paramValues, stabilityValues, 'LineWidth', 2, 'Color', [0.2 0.6 0.8], ...
    'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', [0.2 0.6 0.8]);

% Add labels and grid
title(stabAxes, 'System Stability vs. Parameter Value', 'FontSize', 14, 'FontWeight', 'bold');
xlabel(stabAxes, 'Parameter Value', 'FontSize', 12);
ylabel(stabAxes, 'Stability (1 = Stable, 0 = Unstable)', 'FontSize', 12);
grid(stabAxes, 'on');
ylim(stabAxes, [-0.1 1.1]);

% Add stability transition markers
transitions = findStabilityTransitions(batchResults);
if ~isempty(transitions)
    hold(stabAxes, 'on');
    for t = transitions
        transValue = (paramValues(t) + paramValues(t+1))/2;
        xline(stabAxes, transValue, 'r--', 'LineWidth', 2);
        text(stabAxes, transValue, 0.5, 'Stability Transition', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontWeight', 'bold', 'Color', 'r');
    end
    hold(stabAxes, 'off');
end

% Create stability transition table
if ~isempty(transitions)
    % Create HTML table
    htmlContent = ['<html><head><style>', ...
        'table { border-collapse: collapse; width: 100%; font-family: Arial, sans-serif; }', ...
        'th { background-color: #4472C4; color: white; padding: 8px; text-align: center; font-size: 14px; }', ...
        'td { border: 1px solid #DDDDDD; padding: 8px; text-align: center; font-size: 13px; }', ...
        'tr:nth-child(even) { background-color: #F2F2F2; }', ...
        'tr:hover { background-color: #E8F4F8; }', ...
        '</style></head><body>', ...
        '<table>', ...
        '<tr><th>Transition #</th><th>Parameter Value</th><th>Transition Type</th></tr>'];
    
    for i = 1:length(transitions)
        t = transitions(i);
        % Replace ternary operator with if-else
        if batchResults.stability(t)
            from = 'Stable → Unstable';
        else
            from = 'Unstable → Stable';
        end
        
        % Alternate row color for better readability
        if mod(i, 2) == 0
            bgColor = '#F2F2F2';
        else
            bgColor = '#FFFFFF';
        end
        
        htmlContent = [htmlContent, ...
            sprintf('<tr style="background-color: %s;">', bgColor), ...
            sprintf('<td>%d</td>', i), ...
            sprintf('<td>%.4f</td>', (paramValues(t) + paramValues(t+1))/2), ...
            sprintf('<td>%s</td>', from), ...
            '</tr>'];
    end
    
    htmlContent = [htmlContent, '</table></body></html>'];
    
    % Create HTML component
    uihtml(tab, 'HTMLSource', htmlContent, 'Position', [50 50 900 250]);
else
    % No transitions, display a message
    if all(batchResults.stability)
        msg = 'System remains stable across the entire parameter range.';
        color = [0 0.6 0]; % Green
    else
        msg = 'System is unstable across the entire parameter range.';
        color = [0.8 0 0]; % Red
    end
    
    msgLabel = uilabel(tab, 'Position', [50 150 900 50], 'Text', msg, ...
        'FontSize', 16, 'FontWeight', 'bold', 'FontColor', color, ...
        'HorizontalAlignment', 'center');
end
end