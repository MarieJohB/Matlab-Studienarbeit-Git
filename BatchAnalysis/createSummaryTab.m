function createSummaryTab(tab, batchResults)
% Create a panel for parameter info
infoPanel = uipanel(tab, 'Title', 'Batch Analysis Information', ...
    'Position', [20 400 940 230], 'FontWeight', 'bold');

% Extract parameter info
paramInfo = batchResults.paramInfo;

% Create parameter info text
paramStr = sprintf('%s: %s [%d]', paramInfo.type, paramInfo.coeffType, paramInfo.index);
minVal = paramInfo.min;
maxVal = paramInfo.max;
stepSize = paramInfo.step;
numValues = length(batchResults.paramValues);

% Display parameter info
uilabel(infoPanel, 'Position', [20 180 900 22], 'Text', ...
    sprintf('Parameter Swept: %s', paramStr), 'FontSize', 14);

uilabel(infoPanel, 'Position', [20 150 900 22], 'Text', ...
    sprintf('Range: %.4f to %.4f with step size %.4f (%d values)', ...
    minVal, maxVal, stepSize, numValues), 'FontSize', 14);

% Display transfer function info
G_text = getTransferFunctionText(batchResults.baseG, 'G');
K_text = getTransferFunctionText(batchResults.baseK, 'K');

uilabel(infoPanel, 'Position', [20 100 900 22], 'Text', ...
    sprintf('Base Plant: %s', G_text), 'FontSize', 14);

uilabel(infoPanel, 'Position', [20 70 900 22], 'Text', ...
    sprintf('Base Controller: %s', K_text), 'FontSize', 14);

% Add a key findings section if stability information is available
if isfield(batchResults, 'stability') || isfield(batchResults, 'margins')
    % Create a panel for key findings
    findingsPanel = uipanel(tab, 'Title', 'Key Findings', ...
        'Position', [20 150 940 230], 'FontWeight', 'bold');
    
    vertPosition = 180;
    
    % Add stability findings if available
    if isfield(batchResults, 'stability')
        % Calculate stability transitions
        transitions = findStabilityTransitions(batchResults);
        
        if isempty(transitions)
            if all(batchResults.stability)
                uilabel(findingsPanel, 'Position', [20 vertPosition 900 22], 'Text', ...
                    'System remains stable across the entire parameter range.', ...
                    'FontSize', 14, 'FontColor', [0 0.6 0]);
            elseif ~any(batchResults.stability)
                uilabel(findingsPanel, 'Position', [20 vertPosition 900 22], 'Text', ...
                    'System is unstable across the entire parameter range.', ...
                    'FontSize', 14, 'FontColor', [0.8 0 0]);
            end
        else
            uilabel(findingsPanel, 'Position', [20 vertPosition 900 22], 'Text', ...
                sprintf('Stability Transitions Found at %d parameter values:', length(transitions)), ...
                'FontSize', 14);
            
            for i = 1:min(3, length(transitions))
                idx = transitions(i);
                vertPosition = vertPosition - 25;
                uilabel(findingsPanel, 'Position', [40 vertPosition 900 22], 'Text', ...
                    sprintf('Transition at parameter value = %.4f (index %d)', ...
                    batchResults.paramValues(idx), idx), ...
                    'FontSize', 12, 'FontColor', [0.8 0.4 0]);
            end
        end
        
        vertPosition = vertPosition - 30;
    end
    
    % Add margins findings if available
    if isfield(batchResults, 'margins')
        % Find parameter values with adequate margins
        minGainMargin_dB = 6.0;  % dB
        minPhaseMargin_deg = 30.0;  % degrees
        
        % Initialize arrays for margins
        numPoints = length(batchResults.paramValues);
        gainMargin_dB = zeros(1, numPoints);
        phaseMargin_deg = zeros(1, numPoints);
        
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
            else
                gainMargin_dB(i) = NaN;
                phaseMargin_deg(i) = NaN;
            end
        end
        
        % Count values meeting criteria
        goodGmCount = sum(gainMargin_dB >= minGainMargin_dB & ~isnan(gainMargin_dB));
        goodPmCount = sum(phaseMargin_deg >= minPhaseMargin_deg & ~isnan(phaseMargin_deg));
        bothGoodCount = sum((gainMargin_dB >= minGainMargin_dB & phaseMargin_deg >= minPhaseMargin_deg) & ...
                            ~isnan(gainMargin_dB) & ~isnan(phaseMargin_deg));
        validCount = sum(~isnan(gainMargin_dB) & ~isnan(phaseMargin_deg));
        
        if validCount > 0
            % Find best overall parameter
            normGm = gainMargin_dB / minGainMargin_dB;
            normPm = phaseMargin_deg / minPhaseMargin_deg;
            combinedNorm = normGm + normPm;
            [~, bestIdx] = max(combinedNorm);
            
            % Add margin findings
            uilabel(findingsPanel, 'Position', [20 vertPosition 900 22], 'Text', ...
                sprintf('Margins Analysis: %.1f%% of parameter values have adequate margins.', ...
                bothGoodCount/validCount*100), ...
                'FontSize', 14);
            
            vertPosition = vertPosition - 25;
            
            if bothGoodCount > 0
                uilabel(findingsPanel, 'Position', [40 vertPosition 900 22], 'Text', ...
                    sprintf('Best parameter value for stability margins: %.4f (GM=%.2f dB, PM=%.2f°)', ...
                    batchResults.paramValues(bestIdx), gainMargin_dB(bestIdx), phaseMargin_deg(bestIdx)), ...
                    'FontSize', 12, 'FontColor', [0 0.5 0.7]);
            else
                uilabel(findingsPanel, 'Position', [40 vertPosition 900 22], 'Text', ...
                    'No parameter values have adequate stability margins.', ...
                    'FontSize', 12, 'FontColor', [0.8 0 0]);
            end
        end
    end
end

% Create a quick visualization section
% Plot parameter vs stability if available
if isfield(batchResults, 'stability')
    quickPlotAxes = uiaxes(tab, 'Position', [20 20 940 110]);
    plot(quickPlotAxes, batchResults.paramValues, double(batchResults.stability), ...
        'LineWidth', 2, 'Marker', '.', 'MarkerSize', 15);
    title(quickPlotAxes, 'Parameter vs. Stability');
    xlabel(quickPlotAxes, 'Parameter Value');
    ylabel(quickPlotAxes, 'Stability (0=Unstable, 1=Stable)');
    ylim(quickPlotAxes, [-0.1 1.1]);
    grid(quickPlotAxes, 'on');
end
end