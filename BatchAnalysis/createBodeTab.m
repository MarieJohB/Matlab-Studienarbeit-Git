function createBodeTab(tab, batchResults)
% CREATEBODETAB - Creates the Bode plot visualization tab without panels
% Optimized for 1920x1080 resolution
%
% This function creates a tab showing Bode plots with input fields instead of sliders
% and improved layout for better analysis with elements directly on the tab
%
% Inputs:
%   tab - Parent UI tab object
%   batchResults - Structure containing batch analysis results

% Extract parameter values
paramValues = batchResults.paramValues;

% Create axes for Bode plots directly on tab for fullscreen
magAxes = uiaxes(tab, 'Position', [50 500 1800 450]);
phaseAxes = uiaxes(tab, 'Position', [50 150 1800 350]);

% Create mode selection button group directly on tab
modePanel = uibuttongroup(tab, 'Title', 'Display Mode', ...
    'Position', [100 60 300 60], 'SelectionChangedFcn', @modeChanged, 'FontSize', 14);
singleMode = uiradiobutton(modePanel, 'Text', 'Single Plot', ...
    'Position', [20 5 160 30], 'Value', true, 'FontSize', 14);
compareMode = uiradiobutton(modePanel, 'Text', 'Compare', ...
    'Position', [180 5 160 30], 'FontSize', 14);

% Create dropdown for primary parameter
uilabel(tab, 'Position', [350 80 160 30], 'Text', 'Primary Parameter:', ...
    'FontSize', 14, 'HorizontalAlignment', 'right');
primaryDropdown = uidropdown(tab, ...
    'Position', [520 80 280 30], ...
    'Items', arrayfun(@(x) sprintf('%.6f', x), paramValues, 'UniformOutput', false), ...
    'Value', sprintf('%.6f', paramValues(1)), ...
    'ValueChangedFcn', @(dd, event) updateBodePlot(), 'FontSize', 14);

% Create dropdown for comparison parameter (initially disabled)
uilabel(tab, 'Position', [850 80 180 30], 'Text', 'Comparison Parameter:', ...
    'FontSize', 14, 'HorizontalAlignment', 'right');
comparisonDropdown = uidropdown(tab, ...
    'Position', [1040 80 280 30], ...
    'Items', arrayfun(@(x) sprintf('%.6f', x), paramValues, 'UniformOutput', false), ...
    'Value', sprintf('%.6f', paramValues(min(2, length(paramValues)))), ...
    'Enable', 'off', ...
    'ValueChangedFcn', @(dd, event) updateBodePlot(), 'FontSize', 14);

% Initial plot
updateBodePlot();

% Function to handle mode changes
function modeChanged(~, ~)
    if compareMode.Value
        comparisonDropdown.Enable = 'on';
        compLabel.Visible = 'on';
        legendExplanation.Visible = 'on';
    else
        comparisonDropdown.Enable = 'off';
        compLabel.Visible = 'off';
        legendExplanation.Visible = 'off';
    end
    
    updateBodePlot();
end

% Function to update Bode plots
function updateBodePlot()
    % Get selected parameters
    primaryValue = str2double(primaryDropdown.Value);
    [~, primaryIdx] = min(abs(paramValues - primaryValue));
    
    % Clear plots
    cla(magAxes);
    cla(phaseAxes);
    
    % Get primary Bode data
    bodeData = batchResults.bode{primaryIdx};
    
    if isempty(bodeData)
        % No valid data for this parameter
        text(magAxes, 0.5, 0.5, 'No valid Bode data for primary parameter', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized', ...
            'FontSize', 16, 'FontWeight', 'bold');
        text(phaseAxes, 0.5, 0.5, 'No valid Bode data for primary parameter', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized', ...
            'FontSize', 16, 'FontWeight', 'bold');
        return;
    end
    
    % Extract data
    mag = bodeData.magnitude;
    phase = bodeData.phase;
    omega = bodeData.omega;
    
    % Plot magnitude with improved styling
    semilogx(magAxes, omega, 20*log10(mag), 'LineWidth', 3, ...
        'DisplayName', sprintf('Parameter = %.6f', primaryValue));
    
    % Set properties
    title(magAxes, 'Bode Magnitude Plot', 'FontSize', 16);
    ylabel(magAxes, 'Magnitude (dB)', 'FontSize', 14);
    grid(magAxes, 'on');
    legend(magAxes, 'show', 'Location', 'southwest', 'FontSize', 12);
    
    % Add 0 dB line
    hold(magAxes, 'on');
    yline(magAxes, 0, 'r--', 'LineWidth', 2, 'DisplayName', '0 dB');
    
    % Plot phase with improved styling
    semilogx(phaseAxes, omega, phase, 'LineWidth', 3, ...
        'DisplayName', sprintf('Parameter = %.6f', primaryValue));
    
    % Set properties
    title(phaseAxes, 'Bode Phase Plot', 'FontSize', 16);
    xlabel(phaseAxes, 'Frequency (rad/s)', 'FontSize', 14);
    ylabel(phaseAxes, 'Phase (degrees)', 'FontSize', 14);
    grid(phaseAxes, 'on');
    legend(phaseAxes, 'show', 'Location', 'southwest', 'FontSize', 12);
    
    % Add -180 degree line
    hold(phaseAxes, 'on');
    yline(phaseAxes, -180, 'r--', 'LineWidth', 2, 'DisplayName', '-180°');
    
    % Add comparison plot if enabled
    if compareMode.Value
        comparisonValue = str2double(comparisonDropdown.Value);
        [~, comparisonIdx] = min(abs(paramValues - comparisonValue));
        
        % Skip if trying to compare with itself
        if comparisonIdx ~= primaryIdx
            compData = batchResults.bode{comparisonIdx};
            
            if ~isempty(compData)
                % Plot comparison magnitude
                semilogx(magAxes, compData.omega, 20*log10(compData.magnitude), 'LineWidth', 2.5, ...
                    'LineStyle', '--', 'Color', [0.8 0.4 0], ...
                    'DisplayName', sprintf('Compare = %.6f', comparisonValue));
                
                % Plot comparison phase
                semilogx(phaseAxes, compData.omega, compData.phase, 'LineWidth', 2.5, ...
                    'LineStyle', '--', 'Color', [0.8 0.4 0], ...
                    'DisplayName', sprintf('Compare = %.6f', comparisonValue));
                
                % Try to compute comparison margins
                try
                    [compGm, compPm, compWcg, compWcp] = calcMargins(compData.magnitude, compData.phase, compData.omega);
                    
                    if ~isnan(compGm) && ~isnan(compWcp)
                        % Draw gain margin for comparison (dashed)
                        xline(magAxes, compWcp, 'g--', 'LineWidth', 2, 'DisplayName', ...
                            sprintf('Comp GM: %.2f dB', 20*log10(compGm)));
                    end
                    
                    if ~isnan(compPm) && ~isnan(compWcg)
                        % Draw phase margin for comparison (dashed)
                        xline(phaseAxes, compWcg, 'b--', 'LineWidth', 2, 'DisplayName', ...
                            sprintf('Comp PM: %.2f°', compPm));
                    end
                catch
                    % Skip comparison margins if they can't be calculated
                end
            end
        end
    end
    
    % Add phase and gain margins for primary parameter
    try
        % Extract margins from current Bode data
        [Gm, Pm, Wcg, Wcp] = calcMargins(mag, phase, omega);
        
        if ~isnan(Gm) && ~isnan(Wcp)
            % Draw gain margin
            xline(magAxes, Wcp, 'g-', 'LineWidth', 2, ...
                'DisplayName', sprintf('Gain Margin: %.2f dB', 20*log10(Gm)));
            
            % Mark phase crossover on phase plot
            plot(phaseAxes, Wcp, -180, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', ...
                'DisplayName', 'Phase Crossover');
            
            % Add text annotation with margin values
            text(magAxes, Wcp*1.1, -10, sprintf('GM = %.2f dB', 20*log10(Gm)), ...
                'Color', 'g', 'FontWeight', 'bold', 'FontSize', 12);
        end
        
        if ~isnan(Pm) && ~isnan(Wcg)
            % Draw phase margin
            xline(phaseAxes, Wcg, 'b-', 'LineWidth', 2, ...
                'DisplayName', sprintf('Phase Margin: %.2f°', Pm));
            
            % Mark gain crossover on magnitude plot
            plot(magAxes, Wcg, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', ...
                'DisplayName', 'Gain Crossover');
            
            % Add text annotation with margin values
            text(phaseAxes, Wcg*1.1, -160, sprintf('PM = %.2f°', Pm), ...
                'Color', 'b', 'FontWeight', 'bold', 'FontSize', 12);
        end
    catch
        % Skip margin indications if they can't be calculated
    end
    
    hold(magAxes, 'off');
    hold(phaseAxes, 'off');
end
end