function createJumpTable(tab, paramValues, steadyStateOutput, steadyStateError)
% Create HTML table with jumpability information

% Define threshold for "jumpable" - steady state error near zero
errorThreshold = 0.01;

% Find parameter values where system is "jumpable"
isJumpable = abs(steadyStateError) < errorThreshold;
jumpableCount = sum(isJumpable);
validCount = sum(~isnan(steadyStateError));

if validCount == 0
    % No valid data points
    return;
end

jumpablePercent = jumpableCount / validCount * 100;

% Find parameter ranges where system is jumpable
jumpableRanges = findJumpableRanges(paramValues, isJumpable);

% Create HTML content
htmlContent = ['<html><head><style>', ...
    'table { border-collapse: collapse; width: 100%; font-family: Arial, sans-serif; }', ...
    'th { background-color: #4472C4; color: white; padding: 8px; text-align: center; font-size: 14px; }', ...
    'td { border: 1px solid #DDDDDD; padding: 8px; font-size: 13px; }', ...
    'tr:nth-child(even) { background-color: #F2F2F2; }', ...
    'tr:hover { background-color: #E8F4F8; }', ...
    '.good { color: green; font-weight: bold; }', ...
    '.bad { color: red; font-weight: bold; }', ...
    '.warning { color: orange; font-weight: bold; }', ...
    '</style></head><body>', ...
    '<table>', ...
    '<tr><th colspan="2">Jump Analysis Summary</th></tr>'];

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
            rangeStr = [rangeStr, ', '];
        end
        rangeStr = [rangeStr, sprintf('%.4f - %.4f', jumpableRanges(i, 1), jumpableRanges(i, 2))];
    end
    htmlContent = [htmlContent, '<td class="good">' rangeStr '</td>'];
end

htmlContent = [htmlContent, '</tr>'];

% Add jumpability recommendation
htmlContent = [htmlContent, ...
    '<tr>', ...
    '<td>Recommended parameter value for jumpability</td>'];

if jumpableCount > 0
    % Find parameter value with minimum error
    [~, minErrorIdx] = min(abs(steadyStateError));
    htmlContent = [htmlContent, ...
        sprintf('<td class="good">%.4f (error = %.6f)</td>', ...
        paramValues(minErrorIdx), steadyStateError(minErrorIdx))];
else
    htmlContent = [htmlContent, ...
        '<td class="bad">No jumpable parameter value found</td>'];
end

htmlContent = [htmlContent, '</tr>'];

htmlContent = [htmlContent, '</table></body></html>'];

% Create HTML component
uihtml(tab, 'HTMLSource', htmlContent, 'Position', [50 10 900 30]);
end