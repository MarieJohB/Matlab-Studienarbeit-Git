function updateNyquistPlot(slider, event, nyquistAxes, comparisonAxes, batchResults, valueLabel)
% Get the current parameter index
idx = round(event.Value);

% Update value label
paramValue = batchResults.paramValues(idx);
valueLabel.Text = sprintf('Value = %.4f', paramValue);

% Clear plots
cla(nyquistAxes);
cla(comparisonAxes);

% Get Nyquist data for current parameter
nyquistData = batchResults.nyquist{idx};

if isempty(nyquistData)
    % No valid data for this parameter
    text(nyquistAxes, 0.5, 0.5, 'No valid Nyquist data', ...
        'HorizontalAlignment', 'center', 'Units', 'normalized');
    return;
end

% Extract data
lReal = nyquistData.real;
lImag = nyquistData.imag;
omega = nyquistData.omega;

% Plot Nyquist curve
hold(nyquistAxes, 'on');

% Draw critical point and unit circle for reference
% Determine plot limits based on data
maxAbs = max(sqrt(lReal.^2 + lImag.^2)) * 1.2;
maxAbs = max(maxAbs, 2); % Ensure minimum range

% Draw stability regions
% For Nyquist, stable region depends on whether system has RHP poles
% For simplicity, we'll just highlight the critical point

% Plot main curve
plot(nyquistAxes, lReal, lImag, 'LineWidth', 2, 'Color', [0 0.447 0.741]);

% Mark critical point
plot(nyquistAxes, -1, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
text(nyquistAxes, -1.1, -0.1, 'Critical Point (-1,0)', 'Color', 'r');

% Mark ω=0 and ω=∞ points
plot(nyquistAxes, lReal(1), lImag(1), 'go', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'g');
text(nyquistAxes, lReal(1)+0.1, lImag(1), 'ω=0');

plot(nyquistAxes, lReal(end), lImag(end), 'mo', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'm');
text(nyquistAxes, lReal(end)+0.1, lImag(end), 'ω=∞');

% Draw coordinate axes
plot(nyquistAxes, [-maxAbs maxAbs], [0 0], 'k--');
plot(nyquistAxes, [0 0], [-maxAbs maxAbs], 'k--');

% Add circle around critical point for reference
th = linspace(0, 2*pi, 100);
xunit = cos(th) - 1;
yunit = sin(th);
plot(nyquistAxes, xunit, yunit, 'k:');

% Set axis limits and properties
axis(nyquistAxes, 'equal');
xlim(nyquistAxes, [-maxAbs maxAbs]);
ylim(nyquistAxes, [-maxAbs maxAbs]);
grid(nyquistAxes, 'on');
title(nyquistAxes, sprintf('Nyquist Plot for Parameter = %.4f', paramValue));
xlabel(nyquistAxes, 'Real Axis');
ylabel(nyquistAxes, 'Imaginary Axis');
hold(nyquistAxes, 'off');

% Plot parameter comparison on lower axes
% This will show how the distance to critical point changes with parameter
hold(comparisonAxes, 'on');

% Calculate the minimum distance to critical point for each parameter value
numParams = length(batchResults.paramValues);
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

% Plot the distances
plot(comparisonAxes, batchResults.paramValues, minDistances, 'b-o', 'LineWidth', 2);

% Highlight current parameter
plot(comparisonAxes, paramValue, minDistances(idx), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

% Add critical line at distance = 1 (stability threshold)
yline(comparisonAxes, 1, 'r--', 'LineWidth', 1.5);
text(comparisonAxes, batchResults.paramValues(1), 1.05, 'Critical Distance = 1', 'Color', 'r');

% Set properties
title(comparisonAxes, 'Minimum Distance to Critical Point (-1,0) vs. Parameter');
xlabel(comparisonAxes, 'Parameter Value');
ylabel(comparisonAxes, 'Minimum Distance');
grid(comparisonAxes, 'on');
hold(comparisonAxes, 'off');
end