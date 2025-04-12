function createBatchConfigUI()
% Create figure for batch configuration
batchFig = uifigure('Name', 'Batch Analysis Configuration', 'Position', [300 150 800 700]);
batchFig.Color = [1 1 1];

% Main title
uilabel(batchFig, 'Position', [50 650 700 30], 'Text', 'Control Loop Batch Analysis', ...
    'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold');

% Create panels for different sections
transferFunctionsPanel = uipanel(batchFig, 'Title', 'Transfer Functions', ...
    'Position', [50 450 700 180], 'FontWeight', 'bold');

sweepPanel = uipanel(batchFig, 'Title', 'Parameter Sweep Configuration', ...
    'Position', [50 250 340 180], 'FontWeight', 'bold');

analysisPanel = uipanel(batchFig, 'Title', 'Analyses to Run', ...
    'Position', [410 250 340 180], 'FontWeight', 'bold');

savePanel = uipanel(batchFig, 'Title', 'Save Configuration', ...
    'Position', [50 150 700 80], 'FontWeight', 'bold');

% Transfer function inputs
uilabel(transferFunctionsPanel, 'Position', [20 130 100 22], 'Text', 'Plant G(s):');
gNumLabel = uilabel(transferFunctionsPanel, 'Position', [130 130 100 22], 'Text', 'Numerator:');
gNumEdit = uitextarea(transferFunctionsPanel, 'Position', [230 130 180 22], 'Value', {'1'});

gDenLabel = uilabel(transferFunctionsPanel, 'Position', [130 100 100 22], 'Text', 'Denominator:');
gDenEdit = uitextarea(transferFunctionsPanel, 'Position', [230 100 180 22], 'Value', {'1 1'});

uilabel(transferFunctionsPanel, 'Position', [20 60 100 22], 'Text', 'Controller K(s):');
kNumLabel = uilabel(transferFunctionsPanel, 'Position', [130 60 100 22], 'Text', 'Numerator:');
kNumEdit = uitextarea(transferFunctionsPanel, 'Position', [230 60 180 22], 'Value', {'1'});

kDenLabel = uilabel(transferFunctionsPanel, 'Position', [130 30 100 22], 'Text', 'Denominator:');
kDenEdit = uitextarea(transferFunctionsPanel, 'Position', [230 30 180 22], 'Value', {'1'});

% Preview TF button
previewButton = uibutton(transferFunctionsPanel, 'Text', 'Preview Transfer Functions', ...
    'Position', [450 80 200 30], 'ButtonPushedFcn', @(btn,event) previewTransferFunctions(gNumEdit, gDenEdit, kNumEdit, kDenEdit));

% Parameter sweep configuration
uilabel(sweepPanel, 'Position', [20 140 150 22], 'Text', 'Parameter to Sweep:');
paramDropdown = uidropdown(sweepPanel, 'Position', [170 140 150 22], ...
    'Items', {'G: Numerator [1]', 'G: Denominator [1]', 'G: Denominator [2]', ...
              'K: Numerator [1]', 'K: Denominator [1]'});

% Min, max, step inputs
uilabel(sweepPanel, 'Position', [20 100 100 22], 'Text', 'Min Value:');
minEdit = uieditfield(sweepPanel, 'numeric', 'Position', [170 100 150 22], 'Value', 0.1);

uilabel(sweepPanel, 'Position', [20 70 100 22], 'Text', 'Max Value:');
maxEdit = uieditfield(sweepPanel, 'numeric', 'Position', [170 70 150 22], 'Value', 10);

uilabel(sweepPanel, 'Position', [20 40 100 22], 'Text', 'Step Size:');
stepEdit = uieditfield(sweepPanel, 'numeric', 'Position', [170 40 150 22], 'Value', 0.5);

% Analysis selection
stabCheckbox = uicheckbox(analysisPanel, 'Text', 'Stability Analysis', 'Position', [20 140 150 22], 'Value', true);
nyquistCheckbox = uicheckbox(analysisPanel, 'Text', 'Nyquist Diagram', 'Position', [20 115 150 22], 'Value', true);
bodeCheckbox = uicheckbox(analysisPanel, 'Text', 'Bode Diagram', 'Position', [20 90 150 22], 'Value', true);
marginsCheckbox = uicheckbox(analysisPanel, 'Text', 'Gain/Phase Margins', 'Position', [20 65 150 22], 'Value', true);
keyParamsCheckbox = uicheckbox(analysisPanel, 'Text', 'Key Parameters', 'Position', [20 40 150 22], 'Value', true);
jumpCheckbox = uicheckbox(analysisPanel, 'Text', 'Jump Analysis', 'Position', [20 15 150 22], 'Value', true);

% Additional options
uilabel(savePanel, 'Position', [20 40 100 22], 'Text', 'Save Results To:');
savePathEdit = uieditfield(savePanel, 'text', 'Position', [130 40 440 22], 'Value', '');
browseButton = uibutton(savePanel, 'Text', 'Browse...', 'Position', [580 40 100 22], ...
    'ButtonPushedFcn', @(btn,event) browseSavePath(savePathEdit));

% Action buttons
startButton = uibutton(batchFig, 'Text', 'Start Batch Analysis', 'Position', [300 60 200 50], ...
    'BackgroundColor', [0.3 0.8 0.3], 'FontSize', 14, 'FontWeight', 'bold', ...
    'ButtonPushedFcn', @(btn,event) startBatchAnalysis(batchFig, gNumEdit, gDenEdit, kNumEdit, kDenEdit, ...
                                     paramDropdown, minEdit, maxEdit, stepEdit, ...
                                     stabCheckbox, nyquistCheckbox, bodeCheckbox, marginsCheckbox, keyParamsCheckbox, jumpCheckbox, ...
                                     savePathEdit));

cancelButton = uibutton(batchFig, 'Text', 'Cancel', 'Position', [150 75 120 30], ...
    'BackgroundColor', [0.9 0.3 0.3], ...
    'ButtonPushedFcn', @(btn,event) close(batchFig));

loadButton = uibutton(batchFig, 'Text', 'Load Results', 'Position', [530 75 120 30], ...
    'BackgroundColor', [0.3 0.6 0.9], ...
    'ButtonPushedFcn', @(btn,event) loadBatchResults());
end

% Helper functions for batch configuration UI
function browseSavePath(savePathEdit)
    [file, path] = uiputfile('*.mat', 'Save Batch Analysis Results');
    if isequal(file, 0) || isequal(path, 0)
        % User cancelled
        return;
    end
    savePath = fullfile(path, file);
    savePathEdit.Value = savePath;
end

function previewTransferFunctions(gNumEdit, gDenEdit, kNumEdit, kDenEdit)
    try
        % Get transfer function coefficients
        gNum = str2array(gNumEdit.Value{1});
        gDen = str2array(gDenEdit.Value{1});
        kNum = str2array(kNumEdit.Value{1});
        kDen = str2array(kDenEdit.Value{1});
        
        % Create transfer functions
        G = tf(gNum, gDen);
        K = tf(kNum, kDen);
        
        % Create a preview dialog
        previewFig = uifigure('Name', 'Transfer Function Preview', 'Position', [500 500 400 300]);
        
        % Display G(s)
        [numG, denG] = tfdata(G, 'v');
        gText = sprintf('G(s) = %s / %s', mat2str(numG), mat2str(denG));
        uilabel(previewFig, 'Position', [20 200 360 30], 'Text', 'Plant:', ...
            'FontSize', 14, 'FontWeight', 'bold');
        uilabel(previewFig, 'Position', [20 170 360 30], 'Text', gText, ...
            'FontSize', 12);
        
        % Display K(s)
        [numK, denK] = tfdata(K, 'v');
        kText = sprintf('K(s) = %s / %s', mat2str(numK), mat2str(denK));
        uilabel(previewFig, 'Position', [20 120 360 30], 'Text', 'Controller:', ...
            'FontSize', 14, 'FontWeight', 'bold');
        uilabel(previewFig, 'Position', [20 90 360 30], 'Text', kText, ...
            'FontSize', 12);
        
        % Calculate open-loop transfer function
        L = G * K;
        [numL, denL] = tfdata(L, 'v');
        lText = sprintf('L(s) = %s / %s', mat2str(numL), mat2str(denL));
        uilabel(previewFig, 'Position', [20 40 360 30], 'Text', 'Open-Loop:', ...
            'FontSize', 14, 'FontWeight', 'bold');
        uilabel(previewFig, 'Position', [20 10 360 30], 'Text', lText, ...
            'FontSize', 12);
        
    catch ME
        % Show error message
        uialert(gNumEdit.Parent.Parent, ['Error in transfer function: ' ME.message], 'Invalid Transfer Function', 'error');
    end
end

function array = str2array(str)
    % Convert a space-separated string to a numeric array
    parts = strsplit(str);
    array = zeros(1, length(parts));
    
    for i = 1:length(parts)
        array(i) = str2double(parts{i});
    end
end

function loadBatchResults()
    % Prompt for file
    [file, path] = uigetfile('*.mat', 'Load Batch Analysis Results');
    if isequal(file, 0) || isequal(path, 0)
        % User cancelled
        return;
    end
    
    % Full path
    fullPath = fullfile(path, file);
    
    % Load results
    try
        data = load(fullPath);
        if isfield(data, 'batchResults')
            % Visualize results
            batchVisualization(data.batchResults, fullPath);
        else
            uialert(gcf, 'The selected file does not contain valid batch analysis results.', 'Invalid File', 'error');
        end
    catch ME
        uialert(gcf, ['Error loading file: ' ME.message], 'Load Error', 'error');
    end
end

function startBatchAnalysis(batchFig, gNumEdit, gDenEdit, kNumEdit, kDenEdit, paramDropdown, minEdit, maxEdit, stepEdit, stabCheckbox, nyquistCheckbox, bodeCheckbox, marginsCheckbox, keyParamsCheckbox, jumpCheckbox, savePathEdit)
    % Get all parameter values
    gNum = str2array(gNumEdit.Value{1});
    gDen = str2array(gDenEdit.Value{1});
    kNum = str2array(kNumEdit.Value{1});
    kDen = str2array(kDenEdit.Value{1});
    
    % Create transfer functions
    G = tf(gNum, gDen);
    K = tf(kNum, kDen);
    
    % Get parameter sweep info
    paramStr = paramDropdown.Value;
    [tfType, coeffType, index] = parseParamString(paramStr);
    
    % Get sweep range
    minVal = minEdit.Value;
    maxVal = maxEdit.Value;
    stepSize = stepEdit.Value;
    
    % Create parameter info struct
    paramInfo = struct('type', tfType, 'coeffType', coeffType, 'index', index, ...
                      'min', minVal, 'max', maxVal, 'step', stepSize);
    
    % Get analysis options
    analysisOptions = struct('stability', stabCheckbox.Value, ...
                            'nyquist', nyquistCheckbox.Value, ...
                            'bode', bodeCheckbox.Value, ...
                            'margins', marginsCheckbox.Value, ...
                            'keyParams', keyParamsCheckbox.Value, ...
                            'jump', jumpCheckbox.Value);
    
    % Ensure bode is enabled if margins are requested
    if analysisOptions.margins && ~analysisOptions.bode
        uialert(batchFig, 'Enabling Bode analysis because margins calculation requires it', 'Information', 'info');
        analysisOptions.bode = true;
    end
    
    % Get save path
    savePath = savePathEdit.Value;
    if isempty(savePath)
        % If no path specified, prompt user
        [file, path] = uiputfile('*.mat', 'Save Batch Analysis Results');
        if isequal(file, 0) || isequal(path, 0)
            % User cancelled
            return;
        end
        savePath = fullfile(path, file);
    end
    
    % Close the batch configuration window
    delete(batchFig);
    
    % Run the batch analysis (without trying to use the app)
    batchResults = runBatchAnalysis([], G, K, paramInfo, analysisOptions, savePath);
    
    % Visualize the results
    batchVisualization(batchResults, savePath);
end

function [tfType, coeffType, index] = parseParamString(paramStr)
    % Parse the parameter dropdown selection
    % Example: 'G: Numerator [1]'
    parts = strsplit(paramStr, ' ');
    
    % Get transfer function type (G or K)
    tfType = parts{1}(1); % First character
    
    % Get coefficient type (num or den)
    if contains(lower(parts{2}), 'num')
        coeffType = 'num';
    else
        coeffType = 'den';
    end
    
    % Get index
    indexStr = parts{3}(2:end-1); % Remove brackets
    index = str2double(indexStr);
end