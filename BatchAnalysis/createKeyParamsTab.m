function createKeyParamsTab(tab, batchResults)
% Create key parameter plots across the parameter range

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

% Create 2x2 layout of axes
riseAxes = uiaxes(tab, 'Position', [50 350 425 280]);
settlingAxes = uiaxes(tab, 'Position', [525 350 425 280]);
overshootAxes = uiaxes(tab, 'Position', [50 50 425 280]);
peakAxes = uiaxes(tab, 'Position', [525 50 425 280]);

% Plot rise time
plot(riseAxes, paramValues, riseTime, 'b-o', 'LineWidth', 2);
title(riseAxes, 'Rise Time vs Parameter Value');
xlabel(riseAxes, 'Parameter Value');
ylabel(riseAxes, 'Rise Time (s)');
grid(riseAxes, 'on');

% Plot settling time
plot(settlingAxes, paramValues, settlingTime, 'r-o', 'LineWidth', 2);
title(settlingAxes, 'Settling Time vs Parameter Value');
xlabel(settlingAxes, 'Parameter Value');
ylabel(settlingAxes, 'Settling Time (s)');
grid(settlingAxes, 'on');

% Plot overshoot
plot(overshootAxes, paramValues, overshoot, 'g-o', 'LineWidth', 2);
title(overshootAxes, 'Overshoot vs Parameter Value');
xlabel(overshootAxes, 'Parameter Value');
ylabel(overshootAxes, 'Overshoot (%)');
grid(overshootAxes, 'on');

% Plot peak time
plot(peakAxes, paramValues, peakTime, 'm-o', 'LineWidth', 2);
title(peakAxes, 'Peak Time vs Parameter Value');
xlabel(peakAxes, 'Parameter Value');
ylabel(peakAxes, 'Peak Time (s)');
grid(peakAxes, 'on');

% Add HTML table with threshold indications
createKeyParamsTable(tab, paramValues, riseTime, settlingTime, overshoot, peakTime);
end