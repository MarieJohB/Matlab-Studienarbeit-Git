function updateBodePlot(slider, event, magAxes, phaseAxes, batchResults, valueLabel, compSlider, compValueLabel, compareMode)
% UPDATEBODEPLOT - Updates Bode plots with current parameter value and optional comparison
%
% This function updates the magnitude and phase plots for Bode diagrams,
% with support for comparing two different parameter values
%
% Inputs:
%   slider - Parameter slider UI element
%   event - Slider event data
%   magAxes - Axes for magnitude plot
%   phaseAxes - Axes for phase plot
%   batchResults - Structure containing batch analysis results
%   valueLabel - Label showing current parameter value
%   compSlider - Comparison parameter slider
%   compValueLabel - Label showing comparison parameter value
%   compareMode - Boolean indicating if comparison mode is active

% Get the current parameter index
idx = round(event.Value);

% Update value label
paramValue = batchResults.paramValues(idx);
valueLabel.Text = sprintf('Value = %.4f', paramValue);

% Clear plots
cla(magAxes);
cla(phaseAxes);

% Get Bode data for current parameter
bodeData = batchResults.bode{idx};

if isempty(bodeData)
    % No valid data for this parameter
    text(magAxes, 0.5, 0.5, 'No valid Bode data', ...
        'HorizontalAlignment', 'center', 'Units', 'normalized');
    text(phaseAxes, 0.5, 0.5, 'No valid Bode data', ...
        'HorizontalAlignment', 'center', 'Units', 'normalized');
    return;
end

% Extract data
mag = bodeData.magnitude;
phase = bodeData.phase;
omega = bodeData.omega;

% Plot magnitude
semilogx(magAxes, omega, 20*log10(mag), 'LineWidth', 2, 'DisplayName', sprintf('Parameter = %.4f', paramValue));

% Set properties
title(magAxes, 'Bode Magnitude Plot');
ylabel(magAxes, 'Magnitude (dB)');
grid(magAxes, 'on');
legend(magAxes, 'show', 'Location', 'southwest');

% Add 0 dB line
hold(magAxes, 'on');
yline(magAxes, 0, 'r--', 'LineWidth', 1.5, 'DisplayName', '0 dB');

% Plot phase
semilogx(phaseAxes, omega, phase, 'LineWidth', 2, 'DisplayName', sprintf('Parameter = %.4f', paramValue));

% Set properties
title(phaseAxes, 'Bode Phase Plot');
xlabel(phaseAxes, 'Frequency (rad/s)');
ylabel(phaseAxes, 'Phase (degrees)');
grid(phaseAxes, 'on');
legend(phaseAxes, 'show', 'Location', 'southwest');

% Add -180 degree line
hold(phaseAxes, 'on');
yline(phaseAxes, -180, 'r--', 'LineWidth', 1.5, 'DisplayName', '-180°');

% If in compare mode, add the comparison plot
if compareMode
    compIdx = round(compSlider.Value);
    compValue = batchResults.paramValues(compIdx);
    
    % Skip if trying to compare with itself
    if compIdx ~= idx
        compData = batchResults.bode{compIdx};
        
        if ~isempty(compData)
            % Plot comparison magnitude
            semilogx(magAxes, compData.omega, 20*log10(compData.magnitude), 'LineWidth', 2, ...
                'LineStyle', '--', 'Color', [0.8 0.4 0], ...
                'DisplayName', sprintf('Compare = %.4f', compValue));
            
            % Plot comparison phase
            semilogx(phaseAxes, compData.omega, compData.phase, 'LineWidth', 2, ...
                'LineStyle', '--', 'Color', [0.8 0.4 0], ...
                'DisplayName', sprintf('Compare = %.4f', compValue));
            
            % Draw a little legend in the corner explaining the comparison
            annotation('textbox', [0.65 0.95 0.3 0.04], ...
                'String', 'Solid: Current - Dashed: Comparison', ...
                'EdgeColor', 'none', 'FontSize', 9);
        end
    end
end

% Also plot phase and gain margins if possible
try
    % Extract margins from current Bode data
    [Gm, Pm, Wcg, Wcp] = calcMargins(mag, phase, omega);
    
    if ~isnan(Gm) && ~isnan(Wcp)
        % Draw gain margin
        xline(magAxes, Wcp, 'g-', 'LineWidth', 1.5, 'DisplayName', sprintf('Gain Margin: %.2f dB', 20*log10(Gm)));
        
        % Mark phase crossover on phase plot
        plot(phaseAxes, Wcp, -180, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'Phase Crossover');
    end
    
    if ~isnan(Pm) && ~isnan(Wcg)
        % Draw phase margin
        xline(phaseAxes, Wcg, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('Phase Margin: %.2f°', Pm));
        
        % Mark gain crossover on magnitude plot
        plot(magAxes, Wcg, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'Gain Crossover');
    end
    
    % If in compare mode, also show margins for comparison parameter
    if compareMode && compIdx ~= idx
        compData = batchResults.bode{compIdx};
        if ~isempty(compData)
            try
                [compGm, compPm, compWcg, compWcp] = calcMargins(compData.magnitude, compData.phase, compData.omega);
                
                if ~isnan(compGm) && ~isnan(compWcp)
                    % Draw gain margin for comparison (dashed)
                    xline(magAxes, compWcp, 'g--', 'LineWidth', 1.5, 'DisplayName', ...
                        sprintf('Comp GM: %.2f dB', 20*log10(compGm)));
                end
                
                if ~isnan(compPm) && ~isnan(compWcg)
                    % Draw phase margin for comparison (dashed)
                    xline(phaseAxes, compWcg, 'b--', 'LineWidth', 1.5, 'DisplayName', ...
                        sprintf('Comp PM: %.2f°', compPm));
                end
            catch
                % Skip comparison margins if they can't be calculated
            end
        end
    end
catch
    % Skip margin indications if they can't be calculated
end

hold(magAxes, 'off');
hold(phaseAxes, 'off');
end