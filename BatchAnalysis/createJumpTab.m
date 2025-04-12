function createJumpTab(tab, batchResults)
% Create jump analysis plots across the parameter range

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

% Create 2x1 layout of axes
outputAxes = uiaxes(tab, 'Position', [50 350 900 280]);
errorAxes = uiaxes(tab, 'Position', [50 50 900 280]);

% Plot steady state output
plot(outputAxes, paramValues, steadyStateOutput, 'b-o', 'LineWidth', 2);
title(outputAxes, 'Steady State Output vs Parameter Value');
xlabel(outputAxes, 'Parameter Value');
ylabel(outputAxes, 'Steady State Output');
grid(outputAxes, 'on');

% Add ideal line at y=1 (for unit step)
hold(outputAxes, 'on');
yline(outputAxes, 1, 'r--', 'LineWidth', 1.5);
text(outputAxes, paramValues(1), 1.05, 'Ideal Response (y=1)', 'Color', 'r');
hold(outputAxes, 'off');

% Plot steady state error
plot(errorAxes, paramValues, steadyStateError, 'r-o', 'LineWidth', 2);
title(errorAxes, 'Steady State Error vs Parameter Value');
xlabel(errorAxes, 'Parameter Value');
ylabel(errorAxes, 'Steady State Error');
grid(errorAxes, 'on');

% Add zero error line
hold(errorAxes, 'on');
yline(errorAxes, 0, 'g--', 'LineWidth', 1.5);
text(errorAxes, paramValues(1), 0.05, 'Zero Error', 'Color', 'g');
hold(errorAxes, 'off');

% Create a table with jumpability information
createJumpTable(tab, paramValues, steadyStateOutput, steadyStateError);
end