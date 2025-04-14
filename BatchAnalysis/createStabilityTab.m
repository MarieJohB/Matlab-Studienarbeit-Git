function createStabilityTab(tab, batchResults)
% CREATESTABILITYTAB - Creates the stability visualization tab without panels
% Optimized for 1920x1080 resolution
%
% This function creates a tab showing stability analysis results with a larger
% plot and better organized UI elements directly on the tab

% Get parameter values
paramValues = batchResults.paramValues;
stabilityValues = double(batchResults.stability);


% Create axes for stability plot directly on tab (larger for fullscreen)
stabAxes = uiaxes(tab, 'Position', [50 450 1800 500]);

% Plot stability with improved styling
plot(stabAxes, paramValues, stabilityValues, 'LineWidth', 3, 'Color', [0.2 0.6 0.8], ...
    'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', [0.2 0.6 0.8]);

% Add labels and grid
title(stabAxes, 'System Stability vs. Parameter Value', 'FontSize', 18, 'FontWeight', 'bold');
xlabel(stabAxes, 'Parameter Value', 'FontSize', 14);
ylabel(stabAxes, 'Stability (1 = Stable, 0 = Unstable)', 'FontSize', 14);
grid(stabAxes, 'on');
ylim(stabAxes, [-0.1 1.1]);

% Add stability transition markers with improved visibility
transitions = findStabilityTransitions(batchResults);
if ~isempty(transitions)
    hold(stabAxes, 'on');
    for t = transitions
        transValue = (paramValues(t) + paramValues(t+1))/2;
        xline(stabAxes, transValue, 'r--', 'LineWidth', 2.5);
        text(stabAxes, transValue, 0.5, 'Stability Transition', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontWeight', 'bold', 'FontSize', 14, 'Color', 'r');
    end
    hold(stabAxes, 'off');
end

% Create stability summary section
stableCount = sum(stabilityValues);
totalCount = length(stabilityValues);
stablePercent = stableCount / totalCount * 100;

% Create HTML content for stability info
htmlContent = ['<html><head><style>', ...
    'table { border-collapse: collapse; width: 100%; font-family: Arial, sans-serif; }', ...
    'th { background-color: #4472C4; color: white; padding: 8px; text-align: center; font-size: 16px; }', ...
    'td { border: 1px solid #DDDDDD; padding: 8px; font-size: 14px; }', ...
    'tr:nth-child(even) { background-color: #F2F2F2; }', ...
    'tr:hover { background-color: #E8F4F8; }', ...
    '.header { background-color: #5B9BD5; color: white; font-weight: bold; }', ...
    '.good { color: green; font-weight: bold; }', ...
    '.bad { color: red; font-weight: bold; }', ...
    '.warning { color: orange; font-weight: bold; }', ...
    '</style></head><body>', ...
    '<table>', ...
    '<tr class="header"><th colspan="2">Stability Summary</th></tr>'];

% Add stability percentage
if stablePercent >= 75
    stabilityClass = 'good';
elseif stablePercent >= 25
    stabilityClass = 'warning';
else
    stabilityClass = 'bad';
end

htmlContent = [htmlContent, ...
    '<tr><td>Stable Parameter Values</td><td class="' stabilityClass '">' ...
    sprintf('%.1f%% (%d of %d values)', stablePercent, stableCount, totalCount) ...
    '</td></tr>'];

% Add transitions count
htmlContent = [htmlContent, ...
    '<tr><td>Stability Transitions</td><td>' ...
    sprintf('%d transition(s) found', length(transitions)) ...
    '</td></tr>'];

% Add stability ranges if available
if ~isempty(transitions) || all(stabilityValues == 1) || all(stabilityValues == 0)
    htmlContent = [htmlContent, '<tr><td colspan="2" class="header">Stability Ranges</td></tr>'];
    
    if all(stabilityValues == 1)
        htmlContent = [htmlContent, ...
            '<tr><td colspan="2" class="good">System is stable across the entire parameter range</td></tr>'];
    elseif all(stabilityValues == 0)
        htmlContent = [htmlContent, ...
            '<tr><td colspan="2" class="bad">System is unstable across the entire parameter range</td></tr>'];
    else
        % Find stable ranges
        stableRanges = [];
        currentState = stabilityValues(1);
        rangeStart = paramValues(1);
        
        for i = 2:length(stabilityValues)
            if stabilityValues(i) ~= currentState
                % State changed, record the range
                if currentState == 1
                    stableRanges = [stableRanges; rangeStart, paramValues(i-1)];
                end
                currentState = stabilityValues(i);
                rangeStart = paramValues(i);
            end
        end
        
        % Add the last range if stable
        if currentState == 1
            stableRanges = [stableRanges; rangeStart, paramValues(end)];
        end
        
        if isempty(stableRanges)
            htmlContent = [htmlContent, ...
                '<tr><td colspan="2" class="bad">No stable parameter ranges found</td></tr>'];
        else
            for i = 1:size(stableRanges, 1)
                htmlContent = [htmlContent, ...
                    '<tr><td>Stable Range ' num2str(i) '</td><td class="good">' ...
                    sprintf('%.6f to %.6f', stableRanges(i, 1), stableRanges(i, 2)) ...
                    '</td></tr>'];
            end
        end
    end
end

% Create detailed transition table if there are transitions
if ~isempty(transitions)
    htmlContent = [htmlContent, '<tr><td colspan="2" class="header">Stability Transitions</td></tr>'];
    htmlContent = [htmlContent, '<tr><th>Transition</th><th>Parameter Value</th></tr>'];
    
    for i = 1:length(transitions)
        t = transitions(i);
        transValue = (paramValues(t) + paramValues(t+1))/2;
        
        if stabilityValues(t) == 1
            transType = 'Stable → Unstable';
            transClass = 'bad';
        else
            transType = 'Unstable → Stable';
            transClass = 'good';
        end
        
        htmlContent = [htmlContent, ...
            '<tr><td class="' transClass '">' transType '</td><td>' ...
            sprintf('%.6f (Between %.6f and %.6f)', transValue, paramValues(t), paramValues(t+1)) ...
            '</td></tr>'];
    end
end

htmlContent = [htmlContent, '</table></body></html>'];

% Create HTML component for stability info - directly on tab
stabilityInfoHtml = uihtml(tab, 'HTMLSource', htmlContent, 'Position', [50 50 1800 350]);
end