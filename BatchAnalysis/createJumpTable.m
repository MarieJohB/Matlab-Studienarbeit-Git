function createJumpTable(tab, paramValues, steadyStateOutput, steadyStateError)
    % Create HTML table with zero steady-state error capability information
    
    % Define threshold for zero steady-state error
    errorThreshold = 0.01;
    
    % Find parameter values where system has zero steady-state error
    hasZeroSteadyState = abs(steadyStateError) < errorThreshold;
    zeroSteadyStateCount = sum(hasZeroSteadyState);
    validCount = sum(~isnan(steadyStateError));
    
    if validCount == 0
        % No valid data points
        return;
    end
    
    zeroSteadyStatePercent = zeroSteadyStateCount / validCount * 100;
    
    % Find parameter ranges where system has zero steady-state error
    zeroSteadyStateRanges = findJumpableRanges(paramValues, hasZeroSteadyState);
    
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
        '<tr><th colspan="2">Zero Steady-State Error Analysis Summary</th></tr>'];
    
    % Add zero steady-state percentage
    htmlContent = [htmlContent, ...
        '<tr>', ...
        '<td>Parameter values where system has zero steady-state error (error < ' num2str(errorThreshold) ')</td>'];
    
    if zeroSteadyStatePercent >= 75
        htmlContent = [htmlContent, ...
            sprintf('<td class="good">%.1f%% (%d of %d values)</td>', ...
            zeroSteadyStatePercent, zeroSteadyStateCount, validCount)];
    elseif zeroSteadyStatePercent >= 25
        htmlContent = [htmlContent, ...
            sprintf('<td class="warning">%.1f%% (%d of %d values)</td>', ...
            zeroSteadyStatePercent, zeroSteadyStateCount, validCount)];
    else
        htmlContent = [htmlContent, ...
            sprintf('<td class="bad">%.1f%% (%d of %d values)</td>', ...
            zeroSteadyStatePercent, zeroSteadyStateCount, validCount)];
    end
    
    htmlContent = [htmlContent, '</tr>'];
    
    % Add zero steady-state ranges
    htmlContent = [htmlContent, ...
        '<tr>', ...
        '<td>Parameter ranges where system has zero steady-state error</td>'];
    
    if isempty(zeroSteadyStateRanges)
        htmlContent = [htmlContent, '<td class="bad">No zero steady-state error ranges found</td>'];
    else
        rangeStr = '';
        for i = 1:size(zeroSteadyStateRanges, 1)
            if i > 1
                rangeStr = [rangeStr, ', '];
            end
            rangeStr = [rangeStr, sprintf('%.4f - %.4f', zeroSteadyStateRanges(i, 1), zeroSteadyStateRanges(i, 2))];
        end
        htmlContent = [htmlContent, '<td class="good">' rangeStr '</td>'];
    end
    
    htmlContent = [htmlContent, '</tr>'];
    
    % Add zero steady-state recommendation
    htmlContent = [htmlContent, ...
        '<tr>', ...
        '<td>Recommended parameter value for zero steady-state error</td>'];
    
    if zeroSteadyStateCount > 0
        % Find parameter value with minimum error
        [~, minErrorIdx] = min(abs(steadyStateError));
        htmlContent = [htmlContent, ...
            sprintf('<td class="good">%.4f (error = %.6f)</td>', ...
            paramValues(minErrorIdx), steadyStateError(minErrorIdx))];
    else
        htmlContent = [htmlContent, ...
            '<td class="bad">No parameter value with zero steady-state error found</td>'];
    end
    
    htmlContent = [htmlContent, '</tr>'];
    htmlContent = [htmlContent, '</table></body></html>'];
    
    % Create HTML component
    uihtml(tab, 'HTMLSource', htmlContent, 'Position', [50 10 900 30]);
end