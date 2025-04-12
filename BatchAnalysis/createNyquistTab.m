function createNyquistTab(tab, batchResults)
% Create layout with parameter slider and plot
paramLabel = uilabel(tab, 'Position', [50 630 150 22], 'Text', 'Parameter Value:');

% Create slider for parameter selection
paramValues = batchResults.paramValues;
paramSlider = uislider(tab, 'Position', [200 640 400 3], ...
    'Limits', [1 length(paramValues)], 'Value', 1, ...
    'MajorTicks', 1:ceil(length(paramValues)/5):length(paramValues));

% Create value display
valueLabel = uilabel(tab, 'Position', [620 630 150 22], ...
    'Text', sprintf('Value = %.4f', paramValues(1)));

% Create axes for Nyquist plot
nyquistAxes = uiaxes(tab, 'Position', [50 330 900 300]);

% Create comparison axes for lower plot
comparisonAxes = uiaxes(tab, 'Position', [50 30 900 280]);

% Update function for slider
paramSlider.ValueChangingFcn = @(slider, event) updateNyquistPlot(slider, event, ...
    nyquistAxes, comparisonAxes, batchResults, valueLabel);

% Initial plot
updateNyquistPlot(paramSlider, struct('Value', 1), nyquistAxes, comparisonAxes, batchResults, valueLabel);
end