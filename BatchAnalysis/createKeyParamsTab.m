function createKeyParamsTab(tab, batchResults)
% CREATEKEYPARAMSTAB - Creates the key parameters visualization tab
% with plots directly on tab - optimized for 1920x1080 resolution
%
% Inputs:
%   tab - Parent UI tab object
%   batchResults - Structure containing batch analysis results

% Extract parameter values
paramValues = batchResults.paramValues;
numPoints = length(paramValues);

% Initialize arrays for key parameters
riseTime = zeros(1, numPoints);
settlingTime = zeros(1, numPoints);
overshoot = zeros(1, numPoints);
peakTime = zeros(1, numPoints);

% Extract data for each parameter value
for i = 1:numPoints
    if ~isempty(batchResults.keyParams{i})
        info = batchResults.keyParams{i};
        riseTime(i) = info.RiseTime;
        settlingTime(i) = info.SettlingTime;
        overshoot(i) = info.Overshoot;
        peakTime(i) = info.PeakTime;
    else
        riseTime(i) = NaN;
        settlingTime(i) = NaN;
        overshoot(i) = NaN;
        peakTime(i) = NaN;
    end
end

% Add title label for plots
uilabel(tab, 'Position', [50 970 400 20], 'Text', 'Time Response Parameters', ...
    'FontWeight', 'bold', 'FontSize', 16);

% Create 2x2 layout of axes directly on tab
riseAxes = uiaxes(tab, 'Position', [75 710 850 240]);
settlingAxes = uiaxes(tab, 'Position', [945 710 850 240]);
overshootAxes = uiaxes(tab, 'Position', [75 450 850 240]);
peakAxes = uiaxes(tab, 'Position', [945 450 850 240]);

% Plot rise time with improved styling
plot(riseAxes, paramValues, riseTime, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 8);
title(riseAxes, 'Rise Time vs Parameter Value', 'FontSize', 14);
xlabel(riseAxes, 'Parameter Value', 'FontSize', 12);
ylabel(riseAxes, 'Rise Time (s)', 'FontSize', 12);
grid(riseAxes, 'on');

% Plot settling time with improved styling
plot(settlingAxes, paramValues, settlingTime, 'r-o', 'LineWidth', 2.5, 'MarkerSize', 8);
title(settlingAxes, 'Settling Time vs Parameter Value', 'FontSize', 14);
xlabel(settlingAxes, 'Parameter Value', 'FontSize', 12);
ylabel(settlingAxes, 'Settling Time (s)', 'FontSize', 12);
grid(settlingAxes, 'on');

% Plot overshoot with improved styling
plot(overshootAxes, paramValues, overshoot, 'g-o', 'LineWidth', 2.5, 'MarkerSize', 8);
title(overshootAxes, 'Overshoot vs Parameter Value', 'FontSize', 14);
xlabel(overshootAxes, 'Parameter Value', 'FontSize', 12);
ylabel(overshootAxes, 'Overshoot (%)', 'FontSize', 12);
grid(overshootAxes, 'on');

% Plot peak time with improved styling
plot(peakAxes, paramValues, peakTime, 'm-o', 'LineWidth', 2.5, 'MarkerSize', 8);
title(peakAxes, 'Peak Time vs Parameter Value', 'FontSize', 14);
xlabel(peakAxes, 'Parameter Value', 'FontSize', 12);
ylabel(peakAxes, 'Peak Time (s)', 'FontSize', 12);
grid(peakAxes, 'on');

% Define thresholds for good performance
maxRiseTime = 1.5;     % seconds
maxSettlingTime = 5.0; % seconds
maxOvershoot = 20.0;   % percent

% Count how many parameter values meet the criteria
validCount = sum(~isnan(riseTime));
if validCount == 0
    % No valid data points
    uilabel(tab, 'Position', [750 200 300 40], 'Text', 'No valid time-domain data available', ...
        'FontWeight', 'bold', 'FontSize', 16, 'HorizontalAlignment', 'center');
    return;
end

goodRiseCount = sum(riseTime <= maxRiseTime & ~isnan(riseTime));
goodSettlingCount = sum(settlingTime <= maxSettlingTime & ~isnan(settlingTime));
goodOvershootCount = sum(overshoot <= maxOvershoot & ~isnan(overshoot));

% Calculate percentages
goodRisePercent = goodRiseCount / validCount * 100;
goodSettlingPercent = goodSettlingCount / validCount * 100;
goodOvershootPercent = goodOvershootCount / validCount * 100;

% Find best parameter values
[~, bestRiseIdx] = min(riseTime);
[~, bestSettlingIdx] = min(settlingTime);
[~, bestOvershootIdx] = min(overshoot);

% Create a score based on weighted metrics
weights = [0.3, 0.4, 0.3]; % Rise time, settling time, overshoot
validIndices = find(~isnan(riseTime) & ~isnan(settlingTime) & ~isnan(overshoot));
scores = zeros(1, length(validIndices));

if ~isempty(validIndices)
    % Normalize metrics to 0-1 scale (lower is better)
    normRiseTime = riseTime(validIndices) / max(riseTime(validIndices));
    normSettlingTime = settlingTime(validIndices) / max(settlingTime(validIndices));
    normOvershoot = overshoot(validIndices) / max(max(overshoot(validIndices)), 1);
    
    % Calculate weighted score (lower is better)
    scores = weights(1) * normRiseTime + weights(2) * normSettlingTime + weights(3) * normOvershoot;
    
    % Find best overall parameter
    [~, bestIdx] = min(scores);
    bestParamIdx = validIndices(bestIdx);
    bestParam = paramValues(bestParamIdx);
else
    bestParam = NaN;
    bestParamIdx = NaN;
end

% Create HTML table with threshold information - directly on tab
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
    '<tr class="header"><th>Parameter</th><th>Threshold</th><th>% Meeting Threshold</th><th>Best Value</th><th>At Parameter</th></tr>'];

% Add rise time row
if goodRisePercent > 50
    riseClass = 'good';
elseif goodRisePercent > 25
    riseClass = 'warning';
else
    riseClass = 'bad';
end
htmlContent = [htmlContent, ...
    '<tr>', ...
    '<td>Rise Time</td>', ...
    sprintf('<td>≤ %.2f s</td>', maxRiseTime), ...
    sprintf('<td class="%s">%.1f%%</td>', riseClass, goodRisePercent), ...
    sprintf('<td>%.4f s</td>', riseTime(bestRiseIdx)), ...
    sprintf('<td>%.6f</td>', paramValues(bestRiseIdx)), ...
    '</tr>'];

% Add settling time row
if goodSettlingPercent > 50
    settlingClass = 'good';
elseif goodSettlingPercent > 25
    settlingClass = 'warning';
else
    settlingClass = 'bad';
end
htmlContent = [htmlContent, ...
    '<tr>', ...
    '<td>Settling Time</td>', ...
    sprintf('<td>≤ %.2f s</td>', maxSettlingTime), ...
    sprintf('<td class="%s">%.1f%%</td>', settlingClass, goodSettlingPercent), ...
    sprintf('<td>%.4f s</td>', settlingTime(bestSettlingIdx)), ...
    sprintf('<td>%.6f</td>', paramValues(bestSettlingIdx)), ...
    '</tr>'];

% Add overshoot row
if goodOvershootPercent > 50
    overshootClass = 'good';
elseif goodOvershootPercent > 25
    overshootClass = 'warning';
else
    overshootClass = 'bad';
end
htmlContent = [htmlContent, ...
    '<tr>', ...
    '<td>Overshoot</td>', ...
    sprintf('<td>≤ %.2f%%</td>', maxOvershoot), ...
    sprintf('<td class="%s">%.1f%%</td>', overshootClass, goodOvershootPercent), ...
    sprintf('<td>%.4f%%</td>', overshoot(bestOvershootIdx)), ...
    sprintf('<td>%.6f</td>', paramValues(bestOvershootIdx)), ...
    '</tr>'];

% Add overall best parameter row
if ~isnan(bestParam)
    htmlContent = [htmlContent, ...
        '<tr class="header">', ...
        '<td colspan="5">Overall Best Parameter Value (Weighted Metrics)</td></tr>', ...
        '<tr><td colspan="2">Parameter Value</td>', ...
        sprintf('<td colspan="3">%.6f</td></tr>', bestParam), ...
        '<tr><td>Rise Time</td><td>Settling Time</td><td>Overshoot</td><td>Peak Time</td>', ...
        '<td>Performance Score</td></tr>', ...
        sprintf('<tr><td>%.4f s</td><td>%.4f s</td><td>%.4f%%</td><td>%.4f s</td><td>%.2f/10</td></tr>', ...
            riseTime(bestParamIdx), settlingTime(bestParamIdx), overshoot(bestParamIdx), ...
            peakTime(bestParamIdx), 10*(1-scores(bestIdx)))];
else
    htmlContent = [htmlContent, ...
        '<tr class="header">', ...
        '<td colspan="5">No valid parameter with complete data available for overall score</td></tr>'];
end

htmlContent = [htmlContent, '</table></body></html>'];

% Create HTML component directly on tab
thresholdHtml = uihtml(tab, 'HTMLSource', htmlContent, 'Position', [50 50 1800 350]);
end