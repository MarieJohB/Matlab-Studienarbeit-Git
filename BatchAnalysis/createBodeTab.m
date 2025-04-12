function createBodeTab(tab, batchResults)
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

% Create axes for Bode plots
magAxes = uiaxes(tab, 'Position', [50 350 900 280]);
phaseAxes = uiaxes(tab, 'Position', [50 50 900 280]);

% Update function for slider
paramSlider.ValueChangingFcn = @(slider, event) updateBodePlot(slider, event, ...
    magAxes, phaseAxes, batchResults, valueLabel);

% Initial plot
updateBodePlot(paramSlider, struct('Value', 1), magAxes, phaseAxes, batchResults, valueLabel);
end