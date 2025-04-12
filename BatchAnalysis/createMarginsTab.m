function createMarginsTab(tab, batchResults)
% CREATEMARGINSTAB - Creates the margins visualization tab
%
% This function creates a tab displaying gain and phase margins vs parameter value
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

% Create 2x2 layout of axes
gmAxes = uiaxes(tab, 'Position', [50 350 425 280]);
pmAxes = uiaxes(tab, 'Position', [525 350 425 280]);
wcgAxes = uiaxes(tab, 'Position', [50 50 425 280]);
wcpAxes = uiaxes(tab, 'Position', [525 50 425 280]);

% Plot gain margin
plot(gmAxes, paramValues, gainMargin_dB, 'b-o', 'LineWidth', 2);
title(gmAxes, 'Gain Margin vs Parameter Value');
xlabel(gmAxes, 'Parameter Value');
ylabel(gmAxes, 'Gain Margin (dB)');
grid(gmAxes, 'on');

% Add minimum recommended line (6 dB)
hold(gmAxes, 'on');
yline(gmAxes, 6, 'r--', 'LineWidth', 1.5);
text(gmAxes, paramValues(1), 6.5, 'Minimum Recommended (6 dB)', 'Color', 'r');
hold(gmAxes, 'off');

% Plot phase margin
plot(pmAxes, paramValues, phaseMargin_deg, 'r-o', 'LineWidth', 2);
title(pmAxes, 'Phase Margin vs Parameter Value');
xlabel(pmAxes, 'Parameter Value');
ylabel(pmAxes, 'Phase Margin (degrees)');
grid(pmAxes, 'on');

% Add minimum recommended line (30 degrees)
hold(pmAxes, 'on');
yline(pmAxes, 30, 'r--', 'LineWidth', 1.5);
text(pmAxes, paramValues(1), 32, 'Minimum Recommended (30°)', 'Color', 'r');
hold(pmAxes, 'off');

% Plot gain crossover frequency
plot(wcgAxes, paramValues, gainCrossFreq, 'g-o', 'LineWidth', 2);
title(wcgAxes, 'Gain Crossover Frequency vs Parameter Value');
xlabel(wcgAxes, 'Parameter Value');
ylabel(wcgAxes, 'Frequency (rad/s)');
grid(wcgAxes, 'on');

% Plot phase crossover frequency
plot(wcpAxes, paramValues, phaseCrossFreq, 'm-o', 'LineWidth', 2);
title(wcpAxes, 'Phase Crossover Frequency vs Parameter Value');
xlabel(wcpAxes, 'Parameter Value');
ylabel(wcpAxes, 'Frequency (rad/s)');
grid(wcpAxes, 'on');

% Add HTML table with margin recommendations
createMarginsTable(tab, paramValues, gainMargin_dB, phaseMargin_deg);
end