function createMarginsTable(tab, paramValues, gainMargin_dB, phaseMargin_deg)
% CREATEMARGINSSTABLE - Creates HTML table with stability margin information
%
% Inputs:
%   tab - Parent UI tab object
%   paramValues - Vector of parameter values
%   gainMargin_dB - Vector of gain margins in dB
%   phaseMargin_deg - Vector of phase margins in degrees

% Define thresholds for good stability margins
minGainMargin_dB = 6.0;  % dB
minPhaseMargin_deg = 30.0;  % degrees

% Count how many parameter values meet criteria
validCount = sum(~isnan(gainMargin_dB) & ~isnan(phaseMargin_deg));
if validCount == 0
    % No valid data points
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
[~, bestCompromiseIdx] = max(combinedNorm);

% Create HTML content
htmlContent = ['<html><head><style>', ...
    'table { border-collapse: collapse; width: 100%; font-family: Arial, sans-serif; }', ...
    'th { background-color: #4472C4; color: white; padding: 8px; text-align: center; font-size: 14px; }', ...
    'td { border: 1px solid #DDDDDD; padding: 8px; font-size: 13px; }', ...
    'tr:nth-child(even) { background-color: #F2F2F2; }', ...
    'tr:hover { background-color: #E8F4F8; }', ...
    '.good { color: green; font-weight: bold; }', ...
    '.bad { color: red; font-weight: bold; }', ...
    '</style></head><body>', ...
    '<table>', ...
    '<tr><th>Stability Margin</th><th>Threshold</th><th>% Meeting Threshold</th><th>Best Value</th><th>At Parameter</th></tr>'];

% Add gain margin row
if goodGmPercent > 50
    gmClass = 'good';
else
    gmClass = 'bad';
end
htmlContent = [htmlContent, ...
    '<tr>', ...
    '<td>Gain Margin</td>', ...
    sprintf('<td>≥ %.1f dB</td>', minGainMargin_dB), ...
    sprintf('<td class="%s">%.1f%%</td>', gmClass, goodGmPercent), ...
    sprintf('<td>%.2f dB</td>', bestGm), ...
    sprintf('<td>%.4f</td>', paramValues(bestGmIdx)), ...
    '</tr>'];

% Add phase margin row
if goodPmPercent > 50
    pmClass = 'good';
else
    pmClass = 'bad';
end
htmlContent = [htmlContent, ...
    '<tr>', ...
    '<td>Phase Margin</td>', ...
    sprintf('<td>≥ %.1f°</td>', minPhaseMargin_deg), ...
    sprintf('<td class="%s">%.1f%%</td>', pmClass, goodPmPercent), ...
    sprintf('<td>%.2f°</td>', bestPm), ...
    sprintf('<td>%.4f</td>', paramValues(bestPmIdx)), ...
    '</tr>'];

% Add row for parameter values meeting both criteria
if bothGoodPercent > 50
    bothClass = 'good';
else
    bothClass = 'bad';
end
htmlContent = [htmlContent, ...
    '<tr>', ...
    '<td>Both Margins</td>', ...
    '<td>Both thresholds</td>', ...
    sprintf('<td class="%s">%.1f%%</td>', bothClass, bothGoodPercent)];

if bothGoodCount > 0
    htmlContent = [htmlContent, ...
        sprintf('<td>GM=%.2f dB, PM=%.2f°</td>', ...
        gainMargin_dB(bestCompromiseIdx), phaseMargin_deg(bestCompromiseIdx)), ...
        sprintf('<td>%.4f</td>', paramValues(bestCompromiseIdx))];
else
    htmlContent = [htmlContent, ...
        '<td colspan="2" class="bad">No parameter value meets both criteria</td>'];
end
htmlContent = [htmlContent, '</tr>'];

htmlContent = [htmlContent, '</table></body></html>'];

% Create HTML component
uihtml(tab, 'HTMLSource', htmlContent, 'Position', [200 10 600 30]);
end