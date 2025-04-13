function createMarginsTab(tab, batchResults)
% CREATEMARGINSTAB - Creates the margins visualization tab without panels
% Optimized for 1920x1080 resolution
%
% This function creates a tab displaying gain and phase margins vs parameter value
% directly on the tab without panels
%
% Inputs:
%   tab - Parent UI tab object
%   batchResults - Structure containing batch analysis results

% Extract parameter values
paramValues = batchResults.paramValues;
numPoints = length(paramValues);

% Initialize arrays for margins
gainMargin_dB = zeros(1, numPoints);
phaseMargin_deg = zeros(1, numPoints);
gainCrossFreq = zeros(1, numPoints);
phaseCrossFreq = zeros(1, numPoints);

% Extract data for each parameter value
for i = 1:numPoints
    if ~isempty(batchResults.margins{i})
        margins = batchResults.margins{i};
        
        % Convert gain margin to dB
        if ~isnan(margins.gainMargin)
            gainMargin_dB(i) = 20*log10(margins.gainMargin);
        else
            gainMargin_dB(i) = NaN;
        end
        
        phaseMargin_deg(i) = margins.phaseMargin;
        gainCrossFreq(i) = margins.gainCrossoverFreq;
        phaseCrossFreq(i) = margins.phaseCrossoverFreq;
    else
        gainMargin_dB(i) = NaN;
        phaseMargin_deg(i) = NaN;
        gainCrossFreq(i) = NaN;
        phaseCrossFreq(i) = NaN;
    end
end

% Add title label
uilabel(tab, 'Position', [50 950 400 30], 'Text', 'Stability Margins', ...
    'FontWeight', 'bold', 'FontSize', 16);

% Create axes directly on tab for fullscreen
gmAxes = uiaxes(tab, 'Position', [50 600 850 350]);
pmAxes = uiaxes(tab, 'Position', [920 600 850 350]);
wcgAxes = uiaxes(tab, 'Position', [50 400 850 200]);
wcpAxes = uiaxes(tab, 'Position', [920 400 850 200]);

% Plot gain margin with improved styling
plot(gmAxes, paramValues, gainMargin_dB, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 8);
title(gmAxes, 'Gain Margin vs Parameter Value', 'FontSize', 14);
xlabel(gmAxes, 'Parameter Value', 'FontSize', 12);
ylabel(gmAxes, 'Gain Margin (dB)', 'FontSize', 12);
grid(gmAxes, 'on');

% Add minimum recommended line (6 dB) with better visibility
hold(gmAxes, 'on');
yline(gmAxes, 6, 'r--', 'LineWidth', 2);
text(gmAxes, paramValues(1), 6.5, 'Minimum Recommended (6 dB)', 'Color', 'r', 'FontSize', 12);
hold(gmAxes, 'off');

% Plot phase margin with improved styling
plot(pmAxes, paramValues, phaseMargin_deg, 'r-o', 'LineWidth', 2.5, 'MarkerSize', 8);
title(pmAxes, 'Phase Margin vs Parameter Value', 'FontSize', 14);
xlabel(pmAxes, 'Parameter Value', 'FontSize', 12);
ylabel(pmAxes, 'Phase Margin (degrees)', 'FontSize', 12);
grid(pmAxes, 'on');

% Add minimum recommended line (30 degrees) with better visibility
hold(pmAxes, 'on');
yline(pmAxes, 30, 'r--', 'LineWidth', 2);
text(pmAxes, paramValues(1), 32, 'Minimum Recommended (30°)', 'Color', 'r', 'FontSize', 12);
hold(pmAxes, 'off');

% Plot gain crossover frequency with improved styling
plot(wcgAxes, paramValues, gainCrossFreq, 'g-o', 'LineWidth', 2.5, 'MarkerSize', 8);
title(wcgAxes, 'Gain Crossover Frequency vs Parameter Value', 'FontSize', 14);
xlabel(wcgAxes, 'Parameter Value', 'FontSize', 12);
ylabel(wcgAxes, 'Frequency (rad/s)', 'FontSize', 12);
grid(wcgAxes, 'on');

% Plot phase crossover frequency with improved styling
plot(wcpAxes, paramValues, phaseCrossFreq, 'm-o', 'LineWidth', 2.5, 'MarkerSize', 8);
title(wcpAxes, 'Phase Crossover Frequency vs Parameter Value', 'FontSize', 14);
xlabel(wcpAxes, 'Parameter Value', 'FontSize', 12);
ylabel(wcpAxes, 'Frequency (rad/s)', 'FontSize', 12);
grid(wcpAxes, 'on');

% Add title label for analysis section
uilabel(tab, 'Position', [50 370 400 30], 'Text', 'Stability Margin Analysis', ...
    'FontWeight', 'bold', 'FontSize', 16);

% Define thresholds for good stability margins
minGainMargin_dB = 6.0;  % dB
minPhaseMargin_deg = 30.0;  % degrees

% Count how many parameter values meet criteria
validCount = sum(~isnan(gainMargin_dB) & ~isnan(phaseMargin_deg));
if validCount == 0
    % No valid data points
    uilabel(tab, 'Position', [750 200 300 40], 'Text', 'No valid margins data available', ...
        'FontWeight', 'bold', 'FontSize', 16, 'HorizontalAlignment', 'center');
    return;
end

goodGmCount = sum(gainMargin_dB >= minGainMargin_dB & ~isnan(gainMargin_dB));
goodPmCount = sum(phaseMargin_deg >= minPhaseMargin_deg & ~isnan(phaseMargin_deg));
bothGoodCount = sum((gainMargin_dB >= minGainMargin_dB & phaseMargin_deg >= minPhaseMargin_deg) & ...
                    ~isnan(gainMargin_dB) & ~isnan(phaseMargin_deg));

% Calculate percentages
goodGmPercent = goodGmCount / validCount * 100;
goodPmPercent = goodPmCount / validCount * 100;
bothGoodPercent = bothGoodCount / validCount * 100;

% Find best parameter values
[bestGm, bestGmIdx] = max(gainMargin_dB);
[bestPm, bestPmIdx] = max(phaseMargin_deg);

% Find parameter value with best compromise (highest combined normalized margins)
% Normalize each margin by its minimum requirement
normGm = gainMargin_dB / minGainMargin_dB;
normPm = phaseMargin_deg / minPhaseMargin_deg;
combinedNorm = normGm + normPm;

% Handle NaN values in combinedNorm
validIndices = ~isnan(combinedNorm);
if any(validIndices)
    [bestCombinedScore, bestIdx] = max(combinedNorm(validIndices));
    validIndicesArray = find(validIndices);
    bestCompromiseIdx = validIndicesArray(bestIdx);
    hasBestCompromise = true;
else
    hasBestCompromise = false;
    bestCompromiseIdx = NaN;
    bestCombinedScore = NaN;
end

% Create HTML content similar to key params table
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
    '<tr class="header"><th>Stability Margin</th><th>Threshold</th><th>% Meeting Threshold</th><th>Best Value</th><th>At Parameter</th></tr>'];

% Add gain margin row
if goodGmPercent > 50
    gmClass = 'good';
elseif goodGmPercent > 25
    gmClass = 'warning';
else
    gmClass = 'bad';
end
htmlContent = [htmlContent, ...
    '<tr>', ...
    '<td>Gain Margin</td>', ...
    sprintf('<td>≥ %.1f dB</td>', minGainMargin_dB), ...
    sprintf('<td class="%s">%.1f%%</td>', gmClass, goodGmPercent), ...
    sprintf('<td>%.2f dB</td>', bestGm), ...
    sprintf('<td>%.6f</td>', paramValues(bestGmIdx)), ...
    '</tr>'];

% Add phase margin row
if goodPmPercent > 50
    pmClass = 'good';
elseif goodPmPercent > 25
    pmClass = 'warning';
else
    pmClass = 'bad';
end
htmlContent = [htmlContent, ...
    '<tr>', ...
    '<td>Phase Margin</td>', ...
    sprintf('<td>≥ %.1f°</td>', minPhaseMargin_deg), ...
    sprintf('<td class="%s">%.1f%%</td>', pmClass, goodPmPercent), ...
    sprintf('<td>%.2f°</td>', bestPm), ...
    sprintf('<td>%.6f</td>', paramValues(bestPmIdx)), ...
    '</tr>'];

% Add row for parameter values meeting both criteria
if bothGoodPercent > 50
    bothClass = 'good';
elseif bothGoodPercent > 25
    bothClass = 'warning';
else
    bothClass = 'bad';
end
htmlContent = [htmlContent, ...
    '<tr>', ...
    '<td>Both Margins</td>', ...
    '<td>Both thresholds</td>', ...
    sprintf('<td class="%s">%.1f%%</td>', bothClass, bothGoodPercent)];

if bothGoodCount > 0 && hasBestCompromise
    htmlContent = [htmlContent, ...
        sprintf('<td>GM=%.2f dB, PM=%.2f°</td>', ...
        gainMargin_dB(bestCompromiseIdx), phaseMargin_deg(bestCompromiseIdx)), ...
        sprintf('<td>%.6f</td>', paramValues(bestCompromiseIdx))];
else
    htmlContent = [htmlContent, ...
        '<td colspan="2" class="bad">No parameter value meets both criteria</td>'];
end
htmlContent = [htmlContent, '</tr>'];

% Add detail section for optimal values
if hasBestCompromise
    htmlContent = [htmlContent, ...
        '<tr class="header"><th colspan="5">Optimal Parameter Analysis</th></tr>', ...
        '<tr><td>Parameter with Best Overall Margins</td>', ...
        sprintf('<td colspan="4">%.6f (Combined Score: %.2f)</td></tr>', ...
            paramValues(bestCompromiseIdx), bestCombinedScore), ...
        '<tr><td>Gain Margin</td><td>Phase Margin</td><td>Gain Crossover Frequency</td>', ...
        '<td>Phase Crossover Frequency</td><td>Robustness Assessment</td></tr>', ...
        sprintf('<tr><td>%.2f dB</td><td>%.2f°</td><td>%.4f rad/s</td><td>%.4f rad/s</td>', ...
            gainMargin_dB(bestCompromiseIdx), phaseMargin_deg(bestCompromiseIdx), ...
            gainCrossFreq(bestCompromiseIdx), phaseCrossFreq(bestCompromiseIdx))];
    
    % Add robustness assessment
    if gainMargin_dB(bestCompromiseIdx) >= 10 && phaseMargin_deg(bestCompromiseIdx) >= 45
        htmlContent = [htmlContent, '<td class="good">Excellent</td></tr>'];
    elseif gainMargin_dB(bestCompromiseIdx) >= 6 && phaseMargin_deg(bestCompromiseIdx) >= 30
        htmlContent = [htmlContent, '<td class="good">Good</td></tr>'];
    else
        htmlContent = [htmlContent, '<td class="warning">Marginal</td></tr>'];
    end
    
    % Add recommendations
    if gainMargin_dB(bestCompromiseIdx) < minGainMargin_dB || phaseMargin_deg(bestCompromiseIdx) < minPhaseMargin_deg
        htmlContent = [htmlContent, ...
            '<tr><td colspan="5" class="warning">Controller tuning recommended: Increase margins for better robustness</td></tr>'];
    end
end

htmlContent = [htmlContent, '</table></body></html>'];

% Create HTML component directly on tab
analysisHtml = uihtml(tab, 'HTMLSource', htmlContent, 'Position', [50 50 1800 300]);
end