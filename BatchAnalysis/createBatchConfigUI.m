function createBatchConfigUI()
% Create figure for batch configuration
batchFig = uifigure('Name', 'Control Loop Batch Analysis', 'Position', [300 150 720 550]);
batchFig.Color = [1 1 1];

% Define consistent colors for styling (matching design_controller_auto style)
appColors = struct(...
    'background', [0.95 0.95 0.97], ...     % Light gray background
    'panelHeader', [0.2 0.4 0.7], ...       % Blue header
    'panelBg', [0.95 0.95 0.97], ...        % Light gray panel
    'primary', [0.3 0.6 0.9], ...           % Blue buttons
    'confirm', [0.3 0.8 0.3], ...           % Green confirm button
    'cancel', [0.8 0.3 0.3], ...            % Red cancel button
    'text', [0.2 0.2 0.2], ...              % Dark text
    'lightText', [1 1 1]);                  % White text for dark backgrounds

% Main title panel
titlePanel = uipanel(batchFig, 'Position', [10 520 700 30], ...
    'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
titleLabel = uilabel(titlePanel, 'Position', [0 0 700 30], ...
    'Text', 'Control Loop Batch Analysis', ...
    'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold', ...
    'FontColor', appColors.lightText);

% Create panels for different sections
transferFunctionsPanel = uipanel(batchFig, 'Title', 'Transfer Functions', ...
    'Position', [10 330 700 180], 'FontWeight', 'bold', ...
    'BackgroundColor', appColors.panelBg, 'TitlePosition', 'centertop');

sweepPanel = uipanel(batchFig, 'Title', 'Parameter Sweep Configuration', ...
    'Position', [10 260 700 60], 'FontWeight', 'bold', ...
    'BackgroundColor', appColors.panelBg, 'TitlePosition', 'centertop');

analysisPanel = uipanel(batchFig, 'Title', 'Analyses to Run', ...
    'Position', [10 170 700 80], 'FontWeight', 'bold', ...
    'BackgroundColor', appColors.panelBg, 'TitlePosition', 'centertop');

savePanel = uipanel(batchFig, 'Title', 'Save Configuration', ...
    'Position', [10 100 700 60], 'FontWeight', 'bold', ...
    'BackgroundColor', appColors.panelBg, 'TitlePosition', 'centertop');

% Labels for G(s) and K(s) display boxes
uilabel(transferFunctionsPanel, 'Position', [30 110 50 22], 'Text', 'G(s) =', 'FontWeight', 'bold');
uilabel(transferFunctionsPanel, 'Position', [370 110 50 22], 'Text', 'K(s) =', 'FontWeight', 'bold');

% HTML visualizations inside boxes
gHTMLPreview = uihtml(transferFunctionsPanel, 'Position', [110 85 230 75]);
kHTMLPreview = uihtml(transferFunctionsPanel, 'Position', [450 85 230 75]);

% Input fields for numerator and denominator coefficients
% G(s) inputs
gNumLabel = uilabel(transferFunctionsPanel, 'Position', [20 65 80 22], 'Text', 'Numerator:');
gNumEdit = uitextarea(transferFunctionsPanel, 'Position', [110 65 230 22], 'Value', {'1'});

gDenLabel = uilabel(transferFunctionsPanel, 'Position', [20 35 80 22], 'Text', 'Denominator:');
gDenEdit = uitextarea(transferFunctionsPanel, 'Position', [110 35 230 22], 'Value', {'1 1'});

% K(s) inputs
kNumLabel = uilabel(transferFunctionsPanel, 'Position', [360 65 80 22], 'Text', 'Numerator:');
kNumEdit = uitextarea(transferFunctionsPanel, 'Position', [450 65 230 22], 'Value', {'1'});

kDenLabel = uilabel(transferFunctionsPanel, 'Position', [360 35 80 22], 'Text', 'Denominator:');
kDenEdit = uitextarea(transferFunctionsPanel, 'Position', [450 35 230 22], 'Value', {'1'});

% Update button in the center
previewButton = uibutton(transferFunctionsPanel, 'Text', 'Update', ...
    'Position', [300 5 100 25], 'ButtonPushedFcn', @(btn,event) updateTransferFunctionDisplay(), ...
    'BackgroundColor', appColors.primary, 'FontColor', appColors.lightText);

% Parameter sweep configuration
% Parameter dropdown with more space and better layout
uilabel(sweepPanel, 'Position', [10 10 130 20], 'Text', 'Parameter to Sweep:');
paramDropdown = uidropdown(sweepPanel, 'Position', [135 10 140 20], ...
    'Items', {'G: Numerator [1]', 'G: Denominator [1]', 'G: Denominator [2]', ...
              'K: Numerator [1]', 'K: Denominator [1]'});

% Min, max, step inputs with better layout
uilabel(sweepPanel, 'Position', [295 10 60 20], 'Text', 'Min Value:');
minEdit = uieditfield(sweepPanel, 'numeric', 'Position', [355 10 60 20], 'Value', 0.1);

uilabel(sweepPanel, 'Position', [430 10 60 20], 'Text', 'Max Value:');
maxEdit = uieditfield(sweepPanel, 'numeric', 'Position', [495 10 60 20], 'Value', 10);

uilabel(sweepPanel, 'Position', [570 10 60 20], 'Text', 'Step Size:');
stepEdit = uieditfield(sweepPanel, 'numeric', 'Position', [630 10 60 20], 'Value', 0.5);

% Analysis selection - centered and equally distributed checkboxes
checkboxWidth = 170;
checkboxSpacing = 40;
checkboxTop = 33;
checkboxHeight = 20;

% Calculate positions to center the checkboxes in the panel
panelWidth = 700;
totalCheckboxWidth = 3 * checkboxWidth + 2 * checkboxSpacing;
leftMargin = (panelWidth - totalCheckboxWidth) / 2 + 30;

% First row
stabCheckbox = uicheckbox(analysisPanel, 'Text', 'Stability Analysis', ...
    'Position', [leftMargin checkboxTop checkboxWidth checkboxHeight], 'Value', true);
nyquistCheckbox = uicheckbox(analysisPanel, 'Text', 'Nyquist Diagram', ...
    'Position', [leftMargin + checkboxWidth + checkboxSpacing checkboxTop checkboxWidth checkboxHeight], 'Value', true);
bodeCheckbox = uicheckbox(analysisPanel, 'Text', 'Bode Diagram', ...
    'Position', [leftMargin + 2 * (checkboxWidth + checkboxSpacing) checkboxTop checkboxWidth checkboxHeight], 'Value', true);

% Second row
checkboxTop2 = 10;
marginsCheckbox = uicheckbox(analysisPanel, 'Text', 'Gain/Phase Margins', ...
    'Position', [leftMargin checkboxTop2 checkboxWidth checkboxHeight], 'Value', true);
keyParamsCheckbox = uicheckbox(analysisPanel, 'Text', 'Key Parameters', ...
    'Position', [leftMargin + checkboxWidth + checkboxSpacing checkboxTop2 checkboxWidth checkboxHeight], 'Value', true);
jumpCheckbox = uicheckbox(analysisPanel, 'Text', 'Jump Analysis', ...
    'Position', [leftMargin + 2 * (checkboxWidth + checkboxSpacing) checkboxTop2 checkboxWidth checkboxHeight], 'Value', true);

% Save configuration
uilabel(savePanel, 'Position', [20 10 100 20], 'Text', 'Save Results To:');
savePathEdit = uieditfield(savePanel, 'text', 'Position', [130 10 440 20], 'Value', '');
browseButton = uibutton(savePanel, 'Text', 'Browse...', 'Position', [580 10 100 20], ...
    'ButtonPushedFcn', @(btn,event) browseSavePath(savePathEdit), ...
    'BackgroundColor', appColors.primary, 'FontColor', appColors.lightText);

% Action buttons layout
buttonPanelHeight = 80;

buttonPanel = uipanel(batchFig, 'Title', 'Actions', ...
    'Position', [10 10 700 buttonPanelHeight], 'FontWeight', 'bold', ...
    'BackgroundColor', appColors.panelBg, 'TitlePosition', 'centertop');

cancelButton = uibutton(buttonPanel, 'Text', 'Cancel', 'Position', [120 15 120 30], ...
    'BackgroundColor', appColors.cancel, 'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn,event) close(batchFig));

startButton = uibutton(buttonPanel, 'Text', 'Start Batch Analysis', 'Position', [250 10 200 40], ...
    'BackgroundColor', appColors.confirm, 'FontColor', appColors.lightText, 'FontSize', 14, 'FontWeight', 'bold', ...
    'ButtonPushedFcn', @(btn,event) startBatchAnalysis());

loadButton = uibutton(buttonPanel, 'Text', 'Load Results', 'Position', [460 15 120 30], ...
    'BackgroundColor', appColors.primary, 'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn,event) loadBatchResults());

% Initial rendering of transfer functions
updateTransferFunctionDisplay();

% Helper functions

    % Function to browse for save path
    function browseSavePath(savePathEdit)
        [file, path] = uiputfile('*.mat', 'Save Batch Analysis Results');
        if isequal(file, 0) || isequal(path, 0)
            % User cancelled
            return;
        end
        savePath = fullfile(path, file);
        savePathEdit.Value = savePath;
    end

    % Function to convert string to array of coefficients
    function array = str2array(str)
        % Convert a space-separated string to a numeric array
        parts = strsplit(strrep(str, ',', '.'));
        array = zeros(1, length(parts));
        
        for i = 1:length(parts)
            array(i) = str2double(parts{i});
        end
    end

    % Update the HTML display of transfer functions and parameter dropdown
    function updateTransferFunctionDisplay()
        try
            % Get transfer function coefficients
            gNum = str2array(gNumEdit.Value{1});
            gDen = str2array(gDenEdit.Value{1});
            kNum = str2array(kNumEdit.Value{1});
            kDen = str2array(kDenEdit.Value{1});
            
            % Create transfer functions
            G = tf(gNum, gDen);
            K = tf(kNum, kDen);
            
            % Update G(s) HTML preview
            gHTMLPreview.HTMLSource = createTransferFunctionHTML('G(s)', gNum, gDen);
            
            % Update K(s) HTML preview
            kHTMLPreview.HTMLSource = createTransferFunctionHTML('K(s)', kNum, kDen);
            
            % Update parameter dropdown items based on G and K
            updateParameterDropdown(gNum, gDen, kNum, kDen);
            
        catch ME
            % Show error message - FIXED: Use correct uialert syntax
            uialert(batchFig, ['Error in transfer function: ' ME.message], 'Invalid Transfer Function');
        end
    end

    % Create HTML for displaying a transfer function
    function htmlContent = createTransferFunctionHTML(tfName, num, den)
        % Convert the coefficient vectors to formatted polynomial strings
        numStr = polyToHTMLString(num);
        denStr = polyToHTMLString(den);
        
        % Create an HTML fraction with styling
        htmlContent = ['<html><head><style>', ...
            'body { display: flex; justify-content: center; align-items: center; height: 100%; margin: 0; padding: 0; }', ...
            '.fraction { display: inline-block; vertical-align: middle; margin: 0 auto; text-align: center; }', ...
            '.fraction .num { border-bottom: 1.5px solid black; padding: 5px 5px; font-size: 16px; }', ...
            '.fraction .den { padding: 5px 5px; font-size: 16px; }', ...
            '</style></head><body>', ...
            '<div class="fraction">', ...
            '<div class="num">', numStr, '</div>', ...
            '<div class="den">', denStr, '</div>', ...
            '</div>', ...
            '</body></html>'];
    end

    % Convert polynomial coefficients to HTML string
    function polyStr = polyToHTMLString(coeff)
        % Remove any leading zeros to get the correct degree
        idx = find(coeff ~= 0, 1, 'first');
        if isempty(idx)
            polyStr = '0';
            return;
        end
        coeff = coeff(idx:end);
        
        deg = length(coeff) - 1;  % Determine the degree of the polynomial
        terms = {};  % Cell array to store individual terms
        
        % Iterate over all coefficients
        for i = 1:length(coeff)
            coef = coeff(i);
            exp = deg - (i - 1);  % Determine the exponent for the current term
            
            % Skip zero coefficients as they don't affect the expression
            if coef == 0
                continue;
            end
            
            % Handle signs: add " + " if not the first term
            if coef > 0 && ~isempty(terms)
                term = ' + ';
            elseif coef < 0
                term = ' - ';
                coef = abs(coef);  % Use absolute value for display
            else
                term = '';
            end
            
            % Show the coefficient unless it's 1 and not the constant term (exp==0)
            if coef ~= 1 || exp == 0
                term = strcat(term, num2str(coef));
            end
            
            % Add "s" and superscript for exponents greater than 0
            if exp > 1
                term = strcat(term, 's<sup>', num2str(exp), '</sup>');
            elseif exp == 1
                term = strcat(term, 's');
            end
            
            % Add the formatted term to the cell array
            terms{end + 1} = term;
        end
        
        % Join all terms into a single string
        polyStr = strjoin(terms, '');
        
        % If all coefficients are zero, return "0"
        if isempty(polyStr)
            polyStr = '0';
        end
    end

    % Update parameter dropdown based on current transfer functions
    function updateParameterDropdown(gNum, gDen, kNum, kDen)
        % Create items array for the dropdown
        items = {};
        
        % Add G numerator coefficients
        for i = 1:length(gNum)
            items{end+1} = ['G: Numerator [' num2str(i) ']'];
        end
        
        % Add G denominator coefficients
        for i = 1:length(gDen)
            items{end+1} = ['G: Denominator [' num2str(i) ']'];
        end
        
        % Add K numerator coefficients
        for i = 1:length(kNum)
            items{end+1} = ['K: Numerator [' num2str(i) ']'];
        end
        
        % Add K denominator coefficients
        for i = 1:length(kDen)
            items{end+1} = ['K: Denominator [' num2str(i) ']'];
        end
        
        % Update dropdown
        paramDropdown.Items = items;
        
        % Set default selection if available
        if ~isempty(items)
            paramDropdown.Value = items{1};
        end
    end

    % Load batch results function
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
                % FIXED: Corrected uialert syntax - removed 'error' parameter
                uialert(batchFig, 'The selected file does not contain valid batch analysis results.', 'Invalid File');
            end
        catch ME
            % FIXED: Corrected uialert syntax - removed 'error' parameter
            uialert(batchFig, ['Error loading file: ' ME.message], 'Load Error');
        end
    end

    % Start batch analysis function
    function startBatchAnalysis()
        % Get all parameter values
        try
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
                % FIXED: Corrected uialert syntax - removed 'info' parameter
                uialert(batchFig, 'Enabling Bode analysis because margins calculation requires it', 'Information');
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
            
            % Run the batch analysis
            batchResults = runBatchAnalysis([], G, K, paramInfo, analysisOptions, savePath);
            
            % Visualize the results
            batchVisualization(batchResults, savePath);
            
        catch ME
            % FIXED: Corrected uialert syntax - removed 'error' parameter
            uialert(batchFig, ['Error starting batch analysis: ' ME.message], 'Error');
        end
    end

    % Parse parameter string to get type, coefficient type, and index
    function [tfType, coeffType, index] = parseParamString(paramStr)
        % Parse the parameter dropdown selection
        % Example: 'G: Numerator [1]'
        parts = strsplit(paramStr, {' ', '[', ']'});
        
        % Get transfer function type (G or K)
        tfType = parts{1}(1); % First character
        
        % Get coefficient type (num or den)
        if contains(lower(parts{2}), 'num')
            coeffType = 'num';
        else
            coeffType = 'den';
        end
        
        % Get index
        index = str2double(parts{3});
    end
end