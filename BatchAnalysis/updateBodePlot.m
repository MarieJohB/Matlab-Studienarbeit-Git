function updateBodePlot(slider, event, magAxes, phaseAxes, batchResults, valueLabel)
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
semilogx(magAxes, omega, 20*log10(mag), 'LineWidth', 2);

% Set properties
title(magAxes, sprintf('Bode Magnitude Plot for Parameter = %.4f', paramValue));
ylabel(magAxes, 'Magnitude (dB)');
grid(magAxes, 'on');

% Add 0 dB line
hold(magAxes, 'on');
yline(magAxes, 0, 'r--', 'LineWidth', 1.5);
hold(magAxes, 'off');

% Plot phase
semilogx(phaseAxes, omega, phase, 'LineWidth', 2);

% Set properties
title(phaseAxes, 'Bode Phase Plot');
xlabel(phaseAxes, 'Frequency (rad/s)');
ylabel(phaseAxes, 'Phase (degrees)');
grid(phaseAxes, 'on');

% Add -180 degree line
hold(phaseAxes, 'on');
yline(phaseAxes, -180, 'r--', 'LineWidth', 1.5);
hold(phaseAxes, 'off');

% Also plot phase and gain margins if possible
try
    % Extract margins from current Bode data
    [Gm, Pm, Wcg, Wcp] = calcMargins(mag, phase, omega);
    
    if ~isnan(Gm) && ~isnan(Wcg)
        % Draw gain margin
        hold(magAxes, 'on');
        xline(magAxes, Wcg, 'g-', 'LineWidth', 1.5);
        text(magAxes, Wcg, -5, sprintf('Gain Margin: %.2f dB', 20*log10(Gm)), ...
            'Color', 'g', 'FontWeight', 'bold');
        hold(magAxes, 'off');
        
        % Mark phase crossover on phase plot
        hold(phaseAxes, 'on');
        plot(phaseAxes, Wcp, -180, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
        hold(phaseAxes, 'off');
    end
    
    if ~isnan(Pm) && ~isnan(Wcp)
        % Draw phase margin
        hold(phaseAxes, 'on');
        xline(phaseAxes, Wcp, 'b-', 'LineWidth', 1.5);
        text(phaseAxes, Wcp, -160, sprintf('Phase Margin: %.2f°', Pm), ...
            'Color', 'b', 'FontWeight', 'bold');
        hold(phaseAxes, 'off');
        
        % Mark gain crossover on magnitude plot
        hold(magAxes, 'on');
        plot(magAxes, Wcp, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        hold(magAxes, 'off');
    end
catch
    % Skip margin indications if they can't be calculated
end
end