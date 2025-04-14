function createNyquistTab(tab, batchResults)
% CREATENYQUISTTAB - Creates the Nyquist plot visualization tab with elements directly on tab
%
% Inputs:
%   tab - Parent UI tab object
%   batchResults - Structure containing batch analysis results

% CRITICAL: Setze globale Eigenschaft für unsichtbare Figuren
set(0, 'DefaultFigureVisible', 'off');

% Extract parameter values
paramValues = batchResults.paramValues;

% Create axes for Nyquist plot directly on tab
nyquistAxes = uiaxes(tab, 'Position', [50 410 1700 540]);

% Create comparison axes directly on tab - much larger size for better visibility
comparisonAxes = uiaxes(tab, 'Position', [50 50 1800 340]);

% Create parameter selection controls directly on tab
uilabel(tab, 'Position', [650 20 150 22], 'Text', 'Parameter Value:', ...
    'FontSize', 12, 'HorizontalAlignment', 'right');

% Create dropdown for parameter selection
paramDropdown = uidropdown(tab, ...
    'Position', [820 20 300 30], ...
    'Items', arrayfun(@(x) sprintf('%.6f', x), paramValues, 'UniformOutput', false), ...
    'Value', sprintf('%.6f', paramValues(1)), ...
    'ValueChangedFcn', @(dd, ~) updateNyquistPlot(dd)); % Simplified callback

% Create value display to show more info
valueInfoLabel = uilabel(tab, 'Position', [1130 20 80 22], ...
    'Text', sprintf('(1/%d)', length(paramValues)), ...
    'FontSize', 11, 'HorizontalAlignment', 'left');

% Initial plot
updateNyquistPlot(paramDropdown);

% Stelle DefaultFigureVisible zurück (optional, falls andere Teile der App davon betroffen sein könnten)
% set(0, 'DefaultFigureVisible', 'on');

% Function to update Nyquist plot based on selected parameter
function updateNyquistPlot(dropdown)
    % CRITICAL: Stelle sicher, dass alle neuen Figuren unsichtbar sind
    set(0, 'DefaultFigureVisible', 'off');
    
    % Falls Figure 1 bereits existiert, mache sie explizit unsichtbar
    f1 = findobj(0, 'Type', 'figure', 'Number', 1);
    if ~isempty(f1)
        set(f1, 'Visible', 'off');
    end
    
    % Get the selected parameter value and find corresponding index
    selectedValue = str2double(dropdown.Value);
    [~, idx] = min(abs(paramValues - selectedValue));
    
    % Update the value info label
    valueInfoLabel.Text = sprintf('(%d/%d)', idx, length(paramValues));
    
    % Clear plots
    cla(nyquistAxes);
    cla(comparisonAxes);
    
    % Get Nyquist data for selected parameter
    nyquistData = batchResults.nyquist{idx};
    
    if isempty(nyquistData)
        % No valid data for this parameter
        text(nyquistAxes, 0.5, 0.5, 'No valid Nyquist data', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized', ...
            'FontSize', 16, 'FontWeight', 'bold');
        return;
    end
    
    % Extract data
    lReal = nyquistData.real;
    lImag = nyquistData.imag;
    omega = nyquistData.omega;
    
    % CRITICAL: Explizit den richtigen Axes-Handle verwenden
    axes(nyquistAxes);
    
    % Plot Nyquist curve
    hold(nyquistAxes, 'on');
    
    % Determine plot limits based on data
    maxAbs = max(sqrt(lReal.^2 + lImag.^2)) * 1.2;
    maxAbs = max(maxAbs, 2); % Ensure minimum range
    
    % Plot main curve with improved styling
    plot(nyquistAxes, lReal, lImag, 'LineWidth', 3, 'Color', [0 0.447 0.741]);
    
    % Mark critical point
    plot(nyquistAxes, -1, 0, 'r+', 'MarkerSize', 20, 'LineWidth', 4);
    text(nyquistAxes, -1.1, -0.2, 'Critical Point (-1,0)', 'Color', 'r', ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % Mark ω=0 and ω=∞ points with improved visibility
    plot(nyquistAxes, lReal(1), lImag(1), 'go', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'g');
    text(nyquistAxes, lReal(1)+0.1, lImag(1), 'ω=0', 'FontSize', 14);
    
    plot(nyquistAxes, lReal(end), lImag(end), 'mo', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'm');
    text(nyquistAxes, lReal(end)+0.1, lImag(end), 'ω=∞', 'FontSize', 14);
    
    % Draw coordinate axes
    plot(nyquistAxes, [-maxAbs maxAbs], [0 0], 'k--');
    plot(nyquistAxes, [0 0], [-maxAbs maxAbs], 'k--');
    
    % Add circle around critical point for reference
    th = linspace(0, 2*pi, 100);
    xunit = cos(th) - 1;
    yunit = sin(th);
    plot(nyquistAxes, xunit, yunit, 'k:', 'LineWidth', 1.5);
    
    % Set axis limits and properties
    axis(nyquistAxes, 'equal');
    xlim(nyquistAxes, [-maxAbs maxAbs]);
    ylim(nyquistAxes, [-maxAbs maxAbs]);
    grid(nyquistAxes, 'on');
    title(nyquistAxes, sprintf('Nyquist Plot for Parameter = %.6f', selectedValue), 'FontSize', 16);
    xlabel(nyquistAxes, 'Real Axis', 'FontSize', 14);
    ylabel(nyquistAxes, 'Imaginary Axis', 'FontSize', 14);
    
    % Add stability info based on encirclements of -1+0j
    try
        % Calculate minimum distance to critical point
        distances = sqrt((lReal + 1).^2 + lImag.^2);
        minDist = min(distances);
        
        if minDist < 1
            text(nyquistAxes, 0.05, 0.95, sprintf('Min. distance to critical point: %.4f (< 1)', minDist), ...
                'Units', 'normalized', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 14);
        else
            text(nyquistAxes, 0.05, 0.95, sprintf('Min. distance to critical point: %.4f', minDist), ...
                'Units', 'normalized', 'Color', 'g', 'FontWeight', 'bold', 'FontSize', 14);
        end
    catch
        % Skip distance calculation if it fails
    end
    
    hold(nyquistAxes, 'off');
    
    % Update the distance plot
    updateDistancePlot(idx);
    
    % CRITICAL: Alle zwischengespeicherten Figuren unsichtbar machen
    drawnow; % Lasse alle Grafikoperationen abschließen
    allFigs = findall(0, 'Type', 'figure');
    mainFig = ancestor(tab, 'figure');
    
    for i = 1:length(allFigs)
        if isvalid(allFigs(i)) && allFigs(i) ~= mainFig
            set(allFigs(i), 'Visible', 'off');
        end
    end
end

% Function to update the distance plot (larger size with points)
function updateDistancePlot(currentIdx)
    % Clear plot
    cla(comparisonAxes);
    
    % CRITICAL: Explizit den richtigen Axes-Handle verwenden
    axes(comparisonAxes);
    
    % Calculate the minimum distance to critical point for each parameter value
    numParams = length(paramValues);
    minDistances = zeros(1, numParams);
    
    for i = 1:numParams
        if ~isempty(batchResults.nyquist{i})
            data = batchResults.nyquist{i};
            distances = sqrt((data.real + 1).^2 + data.imag.^2);
            minDistances(i) = min(distances);
        else
            minDistances(i) = NaN;
        end
    end
    
    % Plot the distances with improved styling
    hold(comparisonAxes, 'on');
    
    % Plot with points like in stability tab - larger size for better visibility
    plot(comparisonAxes, paramValues, minDistances, 'LineWidth', 3, ...
        'Color', [0.3 0.6 0.9], 'Marker', '.', 'MarkerSize', 20);
    
    % Highlight current parameter
    plot(comparisonAxes, paramValues(currentIdx), minDistances(currentIdx), 'ro', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'r');
    
    % Add critical line at distance = 1 (stability threshold)
    yline(comparisonAxes, 1, 'r--', 'LineWidth', 2);
    text(comparisonAxes, paramValues(1), 1.05, 'Critical Distance = 1', 'Color', 'r', 'FontSize', 14);
    
    % Set properties with larger fonts for better visibility in the larger plot
    ylabel(comparisonAxes, 'Minimum Distance to Critical Point', 'FontSize', 16);
    xlabel(comparisonAxes, 'Parameter Value', 'FontSize', 16);
    grid(comparisonAxes, 'on');
    xlim(comparisonAxes, [min(paramValues) max(paramValues)]);
    title(comparisonAxes, 'Distance to Critical Point (-1,0) vs. Parameter Value', 'FontSize', 18);
    
    % Add more vertical grid lines for better readability in the larger plot
    ax = gca;
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';
    
    hold(comparisonAxes, 'off');
end
end