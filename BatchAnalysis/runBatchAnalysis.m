function batchResults = runBatchAnalysis(app, G, K, paramInfo, analysisOptions, savePath)
% RUNBATCHANALYSIS - Performs batch analysis of a control system by varying a parameter
%
% This function analyzes a control system by sweeping a parameter through a range
% of values and calculating various performance metrics for each value.
%
% Inputs:
%   app          - App instance (can be empty for command-line usage)
%   G            - Plant transfer function
%   K            - Controller transfer function
%   paramInfo    - Structure with parameter information:
%                  .type - 'G' or 'K' (which transfer function to modify)
%                  .coeffType - 'num' or 'den' (which coefficient type to modify)
%                  .index - Index of coefficient to modify
%                  .min - Minimum parameter value
%                  .max - Maximum parameter value
%                  .step - Parameter step size
%   analysisOptions - Structure with boolean fields for each analysis type:
%                    .stability, .nyquist, .bode, .keyParams, .jump, .margins
%   savePath     - Path where results will be saved
%
% Outputs:
%   batchResults - Structure containing analysis results

% Extract parameter info
tfType = paramInfo.type;        % 'G' or 'K'
coeffType = paramInfo.coeffType; % 'num' or 'den'
index = paramInfo.index;
minVal = paramInfo.min;
maxVal = paramInfo.max;
stepSize = paramInfo.step;

% Create parameter values vector
paramValues = minVal:stepSize:maxVal;
numSteps = length(paramValues);

% Initialize results structure
batchResults = struct();
batchResults.paramInfo = paramInfo;
batchResults.paramValues = paramValues;
batchResults.baseG = G;
batchResults.baseK = K;

% Initialize results containers
if analysisOptions.stability
    batchResults.stability = false(1, numSteps);
end

if analysisOptions.nyquist
    batchResults.nyquist = cell(1, numSteps);
end

if analysisOptions.bode
    batchResults.bode = cell(1, numSteps);
end

if analysisOptions.keyParams
    batchResults.keyParams = cell(1, numSteps);
end

if analysisOptions.jump
    batchResults.jump = cell(1, numSteps);
end

% Add new container for gain and phase margins
if isfield(analysisOptions, 'margins') && analysisOptions.margins
    batchResults.margins = cell(1, numSteps);
end

% Create progress UI
progressFig = uifigure('Name', 'Batch Analysis Progress', 'Position', [400 400 600 200]);

% Create a custom progress bar (compatible with older MATLAB versions)
progressBg = uipanel(progressFig, 'Position', [50 100 500 22], 'BackgroundColor', [0.9 0.9 0.9]);
progressBar = uipanel(progressFig, 'Position', [50 100 50 22], 'BackgroundColor', [0.3 0.8 0.3]);

statusLabel = uilabel(progressFig, 'Position', [50 70 500 22], ...
    'Text', 'Starting batch analysis...');
cancelButton = uibutton(progressFig, 'Text', 'Cancel', 'Position', [250 30 100 30], ...
    'ButtonPushedFcn', @(btn,event) cancelAnalysis());

% Flag for cancellation
cancelled = false;

% Define cancel function
function cancelAnalysis()
    cancelled = true;
    close(progressFig);
end

% Loop through parameter values
for i = 1:numSteps
    % Check if cancelled
    if cancelled
        break;
    end
    
    % Update progress bar
    progress = i / numSteps;
    progressBar.Position = [50 100 progress*500 22]; % Update width based on progress
    statusLabel.Text = sprintf('Processing: Parameter = %.4f (%d/%d)', ...
        paramValues(i), i, numSteps);
    drawnow;
    
    % Modify G or K based on parameter type
    try
        if strcmp(tfType, 'G')
            [num, den] = tfdata(G, 'v');
            if strcmp(coeffType, 'num')
                num(index) = paramValues(i);
            else % 'den'
                den(index) = paramValues(i);
            end
            G_modified = tf(num, den);
            K_modified = K;
        else % 'K'
            [num, den] = tfdata(K, 'v');
            if strcmp(coeffType, 'num')
                num(index) = paramValues(i);
            else % 'den'
                den(index) = paramValues(i);
            end
            K_modified = tf(num, den);
            G_modified = G;
        end
        
        % Create transfer functions for analysis
        [T, S, L, GS, KS] = transferfunctions(G_modified, K_modified);
        
        % Run selected analyses
        if analysisOptions.stability
            batchResults.stability(i) = control_loop_stability(G_modified, K_modified);
        end
        
        if analysisOptions.nyquist
            % Capture Nyquist data without plotting
            [lReal, lImag, omega] = nyquistSiso(L);
            batchResults.nyquist{i} = struct('real', lReal, 'imag', lImag, 'omega', omega);
        end
        
        if analysisOptions.bode
            % Capture Bode data without plotting
            [mag, phase, w] = bode(L, logspace(-2, 3, 200));
            batchResults.bode{i} = struct('magnitude', squeeze(mag), 'phase', squeeze(phase), 'omega', w);
            
            % Calculate margins if that option is enabled
            if isfield(analysisOptions, 'margins') && analysisOptions.margins
                try
                    % Calculate margins from Bode data
                    [Gm, Pm, Wcg, Wcp] = calcMargins(squeeze(mag), squeeze(phase), w);
                    batchResults.margins{i} = struct('gainMargin', Gm, 'phaseMargin', Pm, ...
                                                    'gainCrossoverFreq', Wcg, 'phaseCrossoverFreq', Wcp);
                catch ME_margins
                    warning('Error calculating margins at parameter value %.4f: %s', paramValues(i), ME_margins.message);
                    batchResults.margins{i} = struct('gainMargin', NaN, 'phaseMargin', NaN, ...
                                                    'gainCrossoverFreq', NaN, 'phaseCrossoverFreq', NaN);
                end
            end
        end
        
        if analysisOptions.keyParams
            % Calculate step response parameters
            [y, t] = step(T);
            info = stepinfo(y, t);
            batchResults.keyParams{i} = info;
        end
        
        if analysisOptions.jump
            % Calculate jump parameters (stationary values)
            try
                yss = dcgain(T);
                ess = 1 - yss; % Assuming unit step reference
                batchResults.jump{i} = struct('steadyStateOutput', yss, 'steadyStateError', ess);
            catch
                batchResults.jump{i} = struct('steadyStateOutput', NaN, 'steadyStateError', NaN);
            end
        end
        
    catch ME
        % Log errors but continue with next value
        warning('Error at parameter value %.4f: %s', paramValues(i), ME.message);
        
        % Store NaN/empty for this iteration
        if analysisOptions.stability
            batchResults.stability(i) = false;
        end
        
        if analysisOptions.nyquist
            batchResults.nyquist{i} = [];
        end
        
        if analysisOptions.bode
            batchResults.bode{i} = [];
        end
        
        if isfield(analysisOptions, 'margins') && analysisOptions.margins
            batchResults.margins{i} = [];
        end
        
        if analysisOptions.keyParams
            batchResults.keyParams{i} = [];
        end
        
        if analysisOptions.jump
            batchResults.jump{i} = [];
        end
    end
end

% Close progress UI if still open
if isvalid(progressFig)
    close(progressFig);
end

% Save results
save(savePath, 'batchResults');

% Return results
return;
end