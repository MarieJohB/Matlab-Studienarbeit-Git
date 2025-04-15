function createJumpTab(tab, batchResults)
% CREATEJUMPTAB - Creates jump analysis plots directly on the tab
% optimized for 1920x1080 resolution
%
% Inputs:
%   tab - Parent UI tab object
%   batchResults - Structure containing batch analysis results

% Extract parameter values
paramValues = batchResults.paramValues;
numPoints = length(paramValues);

% Initialize arrays for steady state values
steadyStateOutput = zeros(1, numPoints);
steadyStateError = zeros(1, numPoints);

% Extract data for each parameter value
for i = 1:numPoints
    if ~isempty(batchResults.jump{i})
        jump = batchResults.jump{i};
        steadyStateOutput(i) = jump.steadyStateOutput;
        steadyStateError(i) = jump.steadyStateError;
    else
        steadyStateOutput(i) = NaN;
        steadyStateError(i) = NaN;
    end
end

% Add title label for steady state response
uilabel(tab, 'Position', [50 970 400 20], 'Text', 'Steady State Response', ...
    'FontWeight', 'bold', 'FontSize', 16);

% Create axes directly on tab
outputAxes = uiaxes(tab, 'Position', [50 680 1800 270]);
errorAxes = uiaxes(tab, 'Position', [50 420 1800 270]);

% Plot steady state output with improved styling
plot(outputAxes, paramValues, steadyStateOutput, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 8);
title(outputAxes, 'Steady State Output vs Parameter Value', 'FontSize', 14);
xlabel(outputAxes, 'Parameter Value', 'FontSize', 12);
ylabel(outputAxes, 'Steady State Output', 'FontSize', 12);
grid(outputAxes, 'on');

% Add ideal line at y=1 (for unit step)
hold(outputAxes, 'on');
yline(outputAxes, 1, 'r--', 'LineWidth', 2);
text(outputAxes, paramValues(1), 1.05, 'Ideal Response (y=1)', 'Color', 'r', 'FontSize', 12);
hold(outputAxes, 'off');

% Plot steady state error with improved styling
plot(errorAxes, paramValues, steadyStateError, 'r-o', 'LineWidth', 2.5, 'MarkerSize', 8);
title(errorAxes, 'Steady State Error vs Parameter Value', 'FontSize', 14);
xlabel(errorAxes, 'Parameter Value', 'FontSize', 12);
ylabel(errorAxes, 'Steady State Error', 'FontSize', 12);
grid(errorAxes, 'on');

% Add zero error line
hold(errorAxes, 'on');
yline(errorAxes, 0, 'g--', 'LineWidth', 2);
text(errorAxes, paramValues(1), 0.05, 'Zero Error', 'Color', 'g', 'FontSize', 12);
hold(errorAxes, 'off');

% Define threshold for "jumpable" - steady state error near zero
errorThreshold = 0.01;

% Find parameter values where system is "jumpable"
isJumpable = abs(steadyStateError) < errorThreshold;
jumpableCount = sum(isJumpable);
validCount = sum(~isnan(steadyStateError));

if validCount == 0
    % No valid data points
    uilabel(tab, 'Position', [750 200 300 40], 'Text', 'No valid steady-state data available', ...
        'FontWeight', 'bold', 'FontSize', 16, 'HorizontalAlignment', 'center');
    return;
end

jumpablePercent = jumpableCount / validCount * 100;

% Find parameter ranges where system is jumpable
jumpableRanges = findJumpableRanges(paramValues, isJumpable);

% Create HTML content - enhanced for 1920x1080 directly on tab
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
    '<tr class="header"><th colspan="2">Jump Analysis Summary</th></tr>'];

% Add jumpability percentage
htmlContent = [htmlContent, ...
    '<tr>', ...
    '<td>Parameter values where system is jumpable (error < ' num2str(errorThreshold) ')</td>'];

if jumpablePercent >= 75
    htmlContent = [htmlContent, ...
        sprintf('<td class="good">%.1f%% (%d of %d values)</td>', ...
        jumpablePercent, jumpableCount, validCount)];
elseif jumpablePercent >= 25
    htmlContent = [htmlContent, ...
        sprintf('<td class="warning">%.1f%% (%d of %d values)</td>', ...
        jumpablePercent, jumpableCount, validCount)];
else
    htmlContent = [htmlContent, ...
        sprintf('<td class="bad">%.1f%% (%d of %d values)</td>', ...
        jumpablePercent, jumpableCount, validCount)];
end

htmlContent = [htmlContent, '</tr>'];

% Add jumpable ranges
htmlContent = [htmlContent, ...
    '<tr>', ...
    '<td>Parameter ranges where system is jumpable</td>'];

if isempty(jumpableRanges)
    htmlContent = [htmlContent, '<td class="bad">No jumpable ranges found</td>'];
else
    rangeStr = '';
    for i = 1:size(jumpableRanges, 1)
        if i > 1
            rangeStr = [rangeStr, '<br>'];
        end
        rangeStr = [rangeStr, sprintf('Range %d: %.6f - %.6f', i, jumpableRanges(i, 1), jumpableRanges(i, 2))];
    end
    htmlContent = [htmlContent, '<td class="good">' rangeStr '</td>'];
end

htmlContent = [htmlContent, '</tr>'];

% Add error analysis
if jumpableCount > 0
    % Find parameter values with minimum error
    [minError, minErrorIdx] = min(abs(steadyStateError));
    
    % Find parameter with best output (closest to 1)
    [~, bestOutputIdx] = min(abs(steadyStateOutput - 1));
    
    htmlContent = [htmlContent, ...
        '<tr class="header"><th colspan="2">Optimal Parameter Values</th></tr>', ...
        '<tr><td>Minimum Error Parameter</td>', ...
        sprintf('<td class="good">%.6f (error = %.8f)</td>', ...
            paramValues(minErrorIdx), steadyStateError(minErrorIdx)), ...
        '</tr>', ...
        '<tr><td>Best Output Parameter</td>', ...
        sprintf('<td class="good">%.6f (output = %.6f)</td>', ...
            paramValues(bestOutputIdx), steadyStateOutput(bestOutputIdx)), ...
        '</tr>'];
    
    % Add detail table showing top 5 best parameters
    [sortedErrors, sortIdx] = sort(abs(steadyStateError));
    
    htmlContent = [htmlContent, ...
        '<tr class="header"><th colspan="2">Top 5 Parameters with Lowest Error</th></tr>', ...
        '<tr><th>Parameter Value</th><th>Steady State Error</th></tr>'];
    
    for i = 1:min(5, length(sortIdx))
        idx = sortIdx(i);
        htmlContent = [htmlContent, ...
            '<tr>', ...
            sprintf('<td>%.6f</td>', paramValues(idx)), ...
            sprintf('<td>%.8f</td>', steadyStateError(idx)), ...
            '</tr>'];
    end
    
    % Add recommendations section
    htmlContent = [htmlContent, ...
        '<tr class="header"><th colspan="2">Recommendations</th></tr>'];
    
    if minError < 0.001
        htmlContent = [htmlContent, ...
            '<tr><td colspan="2" class="good">System has excellent steady-state accuracy with parameter = ' ...
            sprintf('%.6f', paramValues(minErrorIdx)) '</td></tr>'];
    elseif minError < 0.01
        htmlContent = [htmlContent, ...
            '<tr><td colspan="2" class="good">System has good steady-state accuracy with parameter = ' ...
            sprintf('%.6f', paramValues(minErrorIdx)) '</td></tr>'];
    else
        htmlContent = [htmlContent, ...
            '<tr><td colspan="2" class="warning">Consider adding integral action to improve steady-state accuracy</td></tr>'];
    end
else
    htmlContent = [htmlContent, ...
        '<tr>', ...
        '<td>Recommended parameter value for jumpability</td>', ...
        '<td class="bad">No jumpable parameter value found</td>', ...
        '</tr>', ...
        '<tr class="header"><th colspan="2">Recommendations</th></tr>', ...
        '<tr><td colspan="2" class="warning">Consider adding integral action to the controller to achieve zero steady-state error</td></tr>'];
end

htmlContent = [htmlContent, '</table></body></html>'];

% Create HTML component directly on tab
jumpHTML = uihtml(tab, 'HTMLSource', htmlContent, 'Position', [50 50 1800 350]);
end