function K = get_user_controller(G)
% GET_USER_CONTROLLER - Create a controller with a consistent UI experience
% This function creates a MATLAB app for selecting a controller type and
% entering its parameters. The UI matches the transfer function UI style
% with consistent panels, colors, and layout.
%
% Parameters:
%   G - Optional plant transfer function for auto-tuning
%
% Returns:
%   K - The controller transfer function

    % Initialize K to empty to ensure it always has a value
    K = [];
    
    % Define UI colors to match transfer function UI
    appColors = struct(...
        'background', [0.95 0.95 0.97], ...       % Light gray background
        'panelHeader', [0.2 0.4 0.7], ...         % Blue panel header
        'panelBg', [0.95 0.95 0.97], ...          % Light panel background
        'buttonPrimary', [0.3 0.6 0.9], ...       % Blue buttons
        'buttonConfirm', [0.3 0.8 0.3], ...       % Green confirm button
        'buttonCancel', [0.8 0.3 0.3], ...        % Red cancel button
        'text', [0.2 0.2 0.2], ...                % Dark text
        'lightText', [1 1 1]);                    % White text for dark backgrounds
    
    % Create main figure with improved styling
    mainFig = uifigure('Name', 'Controller Selection', 'Position', [100 100 700 500]); % Standard height
    mainFig.Color = appColors.background;
    
    % Add title panel
    titlePanel = uipanel(mainFig, 'Position', [10 430 680 60], 'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
    titleLabel = uilabel(titlePanel, 'Text', 'Controller Selection', ...
        'Position', [0 0 680 60], 'FontSize', 20, 'FontWeight', 'bold', ...
        'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');
    
    % Add controller selection panel
    controllerPanel = uipanel(mainFig, 'Title', 'Select Controller Type', ...
        'Position', [10 130 680 290], 'TitlePosition', 'centertop', ...
        'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);
    
    % Create buttons for standard controllers in a grid layout
    % CHANGE: Removed the advanced methods from this grid
    controllers = {'P','PI','PD','PID','PT1','PIT1','I2','PIDT1','Custom'};
    grid_width = 3;
    grid_height = 4;  % Enough for 10 buttons (3 columns * 4 rows = 12 max)
    button_width = 160;
    button_height = 60;
    h_spacing = 30;
    v_spacing = 20;
    
    % Calculate total grid width and height
    total_width = grid_width * button_width + (grid_width - 1) * h_spacing;
    total_height = grid_height * button_height + (grid_height - 1) * v_spacing;
    
    % Calculate starting position to center the grid
    start_x = (680 - total_width) / 2;
    start_y = 185; % Change starting height position
    
    % Create buttons in grid layout
    for i = 1:length(controllers)
        row = ceil(i / grid_width);
        col = mod(i-1, grid_width) + 1;
        
        x_pos = start_x + (col-1) * (button_width + h_spacing);
        y_pos = start_y - (row-1) * (button_height + v_spacing);
        
        % Create button with improved styling
        btn = uibutton(controllerPanel, 'push', 'Text', controllers{i}, ...
            'Position', [x_pos, y_pos, button_width, button_height], ...
            'FontSize', 14, 'FontWeight', 'bold', ...
            'BackgroundColor', appColors.buttonPrimary, ...
            'FontColor', appColors.lightText, ...
            'ButtonPushedFcn', @(~,~) openInputDialog(controllers{i}));
    end
    
    % Add auto-tuning button panel
    autoTunePanel = uipanel(mainFig, 'Position', [10 10 680 110], 'Title', 'Advanced Options', ...
        'TitlePosition', 'centertop', 'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);
    
    % Create auto-tuning button
    if nargin >= 1 && ~isempty(G)
        autoTuneBtn = uibutton(autoTunePanel, 'push', 'Text', 'Automatic Controller Design', ...
            'Position', [190, 15, 300, 60], 'FontSize', 16, 'FontWeight', 'bold', ...
            'BackgroundColor', [0.3 0.6 0.9], 'FontColor', 'white', ...
            'ButtonPushedFcn', @(~,~) openAutoTuningDialog());
    else
        % If G is not provided, show disabled button with explanation
        autoTuneBtn = uibutton(autoTunePanel, 'push', 'Text', 'Automatic Controller Design', ...
            'Position', [190, 40, 300, 40], 'FontSize', 14, ...
            'Enable', 'off', ...
            'ButtonPushedFcn', @(~,~) openAutoTuningDialog());
        
        % Add explanation label
        uilabel(autoTunePanel, 'Text', 'Define plant model G(s) first to enable auto-tuning', ...
            'Position', [0, 15, 680, 20], 'FontSize', 12, 'FontColor', [0.7 0.3 0.3], ...
            'HorizontalAlignment', 'center');
    end
    
    % Handle figure close request
    mainFig.CloseRequestFcn = @handleCloseRequest;
    
    % Block execution until a valid controller is confirmed
    uiwait(mainFig);
    
    %% Nested Functions
    
    % Handle close request for main figure
    function handleCloseRequest(~,~)
        % Ensure K has a value (already initialized to [])
        % This will be the return value when user closes without selecting
        uiresume(mainFig);
        delete(mainFig);
    end
    
    % Modal dialog for controller parameter input
    function openInputDialog(ctrlType)
        % Determine field labels based on controller type
        switch ctrlType
            case 'P'
                fieldLabels = {'Proportional Gain (Kp):'};
                defaultVals = {'1'};
                symbolic = {'Kp'};
                symbolPreview = 'Kp';
            case 'PI'
                fieldLabels = {'Proportional Gain (Kp):', 'Integral Gain (Ki):'};
                defaultVals = {'1', '0.5'};
                symbolic = {'Kp', 'Ki'};
                symbolPreview = 'Kp + Ki/s';
            case 'PD'
                fieldLabels = {'Proportional Gain (Kp):', 'Derivative Gain (Kd):'};
                defaultVals = {'1', '0.1'};
                symbolic = {'Kp', 'Kd'};
                symbolPreview = 'Kp + Kd·s';
            case 'PID'
                fieldLabels = {'Proportional Gain (Kp):', 'Integral Gain (Ki):', 'Derivative Gain (Kd):'};
                defaultVals = {'1', '0.5', '0.1'};
                symbolic = {'Kp', 'Ki', 'Kd'};
                symbolPreview = 'Kp + Ki/s + Kd·s';
            case 'PT1'
                fieldLabels = {'Proportional Gain (Kp):', 'Time Constant (T):'};
                defaultVals = {'1', '1'};
                symbolic = {'Kp', 'T'};
                symbolPreview = 'Kp / (Ts + 1)';
            case 'PIT1'
                fieldLabels = {'Proportional Gain (Kp):', 'Integral Gain (Ki):', 'Time Constant (T):'};
                defaultVals = {'1', '0.5', '1'};
                symbolic = {'Kp', 'Ki', 'T'};
                symbolPreview = '(Kp·s + Ki) / (T·s² + s)';
            case 'I2'
                fieldLabels = {'Integral Gain (Ki):'};
                defaultVals = {'1'};
                symbolic = {'Ki'};
                symbolPreview = 'Ki / s²';
            case 'PIDT1'
                fieldLabels = {'Proportional Gain (Kp):', 'Integral Gain (Ki):', 'Derivative Gain (Kd):', 'Time Constant (T):'};
                defaultVals = {'1', '0.5', '0.1', '1'};
                symbolic = {'Kp', 'Ki', 'Kd', 'T'};
                symbolPreview = '(Kd·s² + Kp·s + Ki) / (T·s² + s)';
            case 'Custom'
                fieldLabels = {'Numerator coefficients [b0 b1 ...]:', 'Denominator coefficients [a0 a1 ...]:'};
                defaultVals = {'1 0', '1 1'};
                symbolic = {'num', 'den'};
                symbolPreview = 'numerator / denominator';
            case 'Pole Placement'
                fieldLabels = {'Desired Bandwidth (rad/s):', 'Damping Ratio:', 'Derivative Filter (epsilon):'};
                defaultVals = {'0.3', '1.0', '0.1'}; % Safer default values for unstable systems
                symbolic = {'ω', 'ζ', 'ε'};
                symbolPreview = 'Pole Placement: ω = desired bandwidth, ζ = damping ratio';
        end
        numFields = length(fieldLabels);
        
        % Create a modal UI figure with styling to match transfer function UI
        dlg = uifigure('Name', [ctrlType ' Controller Parameters'], ...
                     'Position', [200 200 600 450], 'WindowStyle', 'modal');
        dlg.Color = appColors.background;
        
        % Title panel
        titlePnl = uipanel(dlg, 'Position', [10 390 580 50], 'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
        titleLbl = uilabel(titlePnl, 'Text', [ctrlType ' Controller Configuration'], ...
            'Position', [0 0 580 50], 'FontSize', 16, 'FontWeight', 'bold', ...
            'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');
        
        % Preview panel
        previewPanel = uipanel(dlg, 'Title', 'Controller Preview', ...
            'Position', [10 280 580 100], 'TitlePosition', 'centertop', ...
            'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);
        
        % Create UIHTML component for the preview with the same styling as transfer function UI
        previewHTML = uihtml(previewPanel, 'Position', [10 5 560 70], ...
            'HTMLSource', getSymbolicPreviewHTML(ctrlType, symbolPreview));
        
        % Parameters panel - adjust height based on number of fields
        paramHeight = 40 * numFields + 20;  % Dynamic height based on number of fields
        paramsPanel = uipanel(dlg, 'Title', 'Controller Parameters', ...
            'Position', [10 270 - paramHeight, 580, paramHeight], 'TitlePosition', 'centertop', ...
            'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);
        
        % Create labels and input fields for controller parameters
        fields = cell(1, numFields);
        for j = 1:numFields
            yPos = paramHeight - (j * 40) - 15;
            
            % Create label with symbolic name shown
            uilabel(paramsPanel, 'Position', [20 yPos 250 30], ...
                'Text', [fieldLabels{j} ' (' symbolic{j} ')'], ...
                'FontSize', 12);
            
            % Create edit field with default value
            fields{j} = uieditfield(paramsPanel, 'text', ...
                'Position', [280 yPos 280 30], ...
                'Value', defaultVals{j}, 'FontSize', 12);
        end
        
        % Buttons panel with no border
        btnPanel = uipanel(dlg, 'Position', [10 10 580 70], 'BackgroundColor', appColors.background);
        
        % Create buttons with same style as in transfer function UI (no individual borders)
        previewBtn = uibutton(btnPanel, 'push', 'Text', 'Preview', ...
            'Position', [120 15 100 40], 'FontSize', 14, ...
            'BackgroundColor', appColors.buttonPrimary, ...
            'FontColor', appColors.lightText, ...
            'ButtonPushedFcn', @(~,~) previewCallbackK());
        
        confirmBtn = uibutton(btnPanel, 'push', 'Text', 'Confirm', ...
            'Position', [240 15 100 40], 'FontSize', 14, ...
            'BackgroundColor', appColors.buttonConfirm, ...
            'FontColor', appColors.lightText, ...
            'ButtonPushedFcn', @(~,~) confirmCallbackK());
        
        cancelBtn = uibutton(btnPanel, 'push', 'Text', 'Cancel', ...
            'Position', [360 15 100 40], 'FontSize', 14, ...
            'BackgroundColor', appColors.buttonCancel, ...
            'FontColor', appColors.lightText, ...
            'ButtonPushedFcn', @(~,~) cancelCallbackK());
            
        % Handle figure close request
        dlg.CloseRequestFcn = @cancelCallbackK;
        
        % Callback for Preview: update the HTML preview area
        function previewCallbackK()
            % If any input field is empty, show the symbolic preview
            if any(cellfun(@(x) isempty(x.Value), fields))
                previewHTML.HTMLSource = getSymbolicPreviewHTML(ctrlType, symbolPreview);
                return;
            end
            
            % Otherwise, try to compute numeric transfer function preview
            try
                % Special handling for Pole Placement
                if strcmp(ctrlType, 'Pole Placement')
                    try
                        % Simple plant analysis for preview
                        p = pole(G);
                        try
                            z = zero(G);
                        catch
                            z = [];
                        end
                        
                        % Create a basic plantInfo structure
                        plantInfo = struct();
                        plantInfo.isUnstable = any(real(p) > 0);
                        plantInfo.poles = p;
                        plantInfo.zeros = z;
                    catch
                        % If analysis fails, assume a stable plant
                        plantInfo = struct();
                        plantInfo.isUnstable = false;
                    end
                    % Extract parameters
                    bandwidth = convertValue(fields{1}.Value);
                    damping = convertValue(fields{2}.Value);
                    epsilon = convertValue(fields{3}.Value);
                    
                    % Create a representative PID controller for preview
                    % This is just for preview purposes and doesn't need the actual pole placement calculation
                    if plantInfo.isUnstable
                        % For unstable plants, use more conservative parameters
                        Kp = 1.0;
                        Ki = 0.01;
                        Kd = 5.0;
                    else
                        % For stable plants
                        Kp = bandwidth;
                        Ki = bandwidth^2 / 10;
                        Kd = bandwidth / 5;
                    end
                    
                    % Create the HTML for the fraction representing a PID controller
                    numHTML = sprintf('%.2fs<sup>2</sup> + %.2fs + %.3f', Kd, Kp, Ki);
                    denHTML = sprintf('%.3fs<sup>2</sup> + s', epsilon*Kd);
                    
                    % Create HTML fraction with exact same styling as other controllers
                    formulaHTML = ['<html><head><style>', ...
                        'body { display: flex; justify-content: center; align-items: center; height: 100%; margin: 0; padding: 0; }', ...
                        '.fraction { display: inline-block; vertical-align: middle; margin: 0 auto; text-align: center; }', ...
                        '.fraction .num { border-bottom: 1px solid black; padding: 5px 15px; font-size: 16px; }', ...
                        '.fraction .den { padding: 5px 15px; font-size: 16px; }', ...
                        '</style></head><body>', ...
                        '<div class="fraction">', ...
                        '<div class="num">', numHTML, '</div>', ...
                        '<div class="den">', denHTML, '</div>', ...
                        '</div>', ...
                        '<div style="margin-top: 10px; font-size: 12px; color: #666;">PID approximation (preview only)</div>', ...
                        '</body></html>'];
                    
                    previewHTML.HTMLSource = formulaHTML;
                    return;
                end
                
                % For other controller types, compute normally
                Ktemp = computeController(ctrlType, fields);
                
                % Extract numerator and denominator from the transfer function (assumed SISO)
                numCell = Ktemp.Numerator;
                denCell = Ktemp.Denominator;
                if isempty(numCell) || isempty(denCell)
                    previewHTML.HTMLSource = '<html><body style="font-size:16px; text-align:center;">Preview: Invalid input.</body></html>';
                    return;
                end
                numVec = numCell{1};
                denVec = denCell{1};
                % Convert polynomial coefficients to HTML-formatted strings without asterisks
                numHTML = polyToHTMLString(numVec);
                denHTML = polyToHTMLString(denVec);
                
                % Create HTML fraction with exact same styling as transfer function UI
                formulaHTML = ['<html><head><style>', ...
                    'body { display: flex; justify-content: center; align-items: center; height: 100%; margin: 0; padding: 0; }', ...
                    '.fraction { display: inline-block; vertical-align: middle; margin: 0 auto; text-align: center; }', ...
                    '.fraction .num { border-bottom: 1px solid black; padding: 5px 15px; font-size: 16px; }', ...
                    '.fraction .den { padding: 5px 15px; font-size: 16px; }', ...
                    '</style></head><body>', ...
                    '<div class="fraction">', ...
                    '<div class="num">', numHTML, '</div>', ...
                    '<div class="den">', denHTML, '</div>', ...
                    '</div>', ...
                    '</body></html>'];
                
                previewHTML.HTMLSource = formulaHTML;
            catch ME
                disp(['Preview error: ' ME.message]);
                previewHTML.HTMLSource = '<html><body style="font-size:16px; text-align:center; color:red;">Preview: Invalid input. ' + ME.message + '</body></html>';
            end
        end
        
        % Callback for Confirm: validate input, create transfer function and close dialog & main figure
        function confirmCallbackK()
            try
                K = computeController(ctrlType, fields);
                % Resume and close safely
                uiresume(dlg);
                delete(dlg);
                uiresume(mainFig);
                delete(mainFig);
            catch ME
                disp(['Error confirming controller: ' ME.message]);
                uialert(dlg, ['Invalid input: ' ME.message], 'Error');
            end
        end
        
        % Callback for Cancel: close the dialog and return to the controller selection window
        function cancelCallbackK()
            uiresume(dlg);
            delete(dlg);
        end
        
        uiwait(dlg);
    end

    % Opens the Auto-Tuning dialog with matching styling
    function openAutoTuningDialog()
    % Use the G plant model passed as parameter to the outer function
    if isempty(G)
         uialert(mainFig, 'Please define plant model G(s) first.', 'Plant Missing');
         return;
    end
    
    % Check stability and order of G
    isStable = all(real(pole(G)) < 0);
    plantOrder = length(pole(G));
    
    % Create a figure with improved styling
    autoTuneFig = uifigure('Name', 'Automatic Controller Design', 'Position', [300 150 700 720]); % Increased height
    autoTuneFig.Color = appColors.background;
    
    % Add a professional-looking title and plant info
    titlePanel = uipanel(autoTuneFig, 'Position', [10 650 680 60], 'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
    titleLabel = uilabel(titlePanel, 'Text', 'Controller Auto-Tuning', ...
        'Position', [0 0 680 60], 'FontSize', 20, 'FontWeight', 'bold', ...
        'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');
    
    % Plant info display with improved styling
    if isStable
        stabilityText = 'stable';
        stabilityColor = [0.2 0.7 0.3]; % Green for stable
    else
        stabilityText = 'unstable';
        stabilityColor = [0.7 0.3 0.2]; % Red for unstable
    end
    
    plantInfoText = sprintf('Plant G(s): %s, Order %d', stabilityText, plantOrder);
    plantInfoLabel = uilabel(autoTuneFig, 'Text', plantInfoText, ...
        'Position', [10 620 680 20], 'FontWeight', 'bold', ...
        'FontColor', stabilityColor, 'HorizontalAlignment', 'center');
    
    % Create main panels with consistent height and spacing
    % All panels now have consistent positioning and proportions
    designPanel = uipanel(autoTuneFig, 'Title', 'Design Methods', ...
        'Position', [10 460 680 155], 'TitlePosition', 'centertop', ...
        'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);
    
    controllerPanel = uipanel(autoTuneFig, 'Title', 'Controller Configuration', ...
        'Position', [10 390 680 60], 'TitlePosition', 'centertop', ...
        'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);
    
    performancePanel = uipanel(autoTuneFig, 'Title', 'Performance Requirements', ...
        'Position', [10 210 680 170], 'TitlePosition', 'centertop', ...
        'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);
        
    resultsPanel = uipanel(autoTuneFig, 'Title', 'Results', ...
        'Position', [10 10 680 190], 'TitlePosition', 'centertop', ...
        'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

    
    % Methods selection with better layout
    methods = {
        'Ziegler-Nichols (Oscillation)', isApplicable(G, 'ZN-Oscillation'),
        'Ziegler-Nichols (Step)', isApplicable(G, 'ZN-Step'),
        'Aström', isApplicable(G, 'Astrom'),
        'Loop-Shaping', isApplicable(G, 'Loop-Shaping'),
        'IMC (Internal Model Control)', isApplicable(G, 'IMC'),
        'MIGO (M-constrained Integral Gain Optimization)', isApplicable(G, 'MIGO'),
        'Pole Placement', isApplicable(G, 'Pole-Placement'),
        'Compensation Controller', isApplicable(G, 'Compensation-Controller')
    };
    
    % Design Methods panel - perfectly centered elements with consistent spacing
    numMethods = size(methods, 1);
    numRows = ceil(numMethods / 2);
    
    % Panel dimensions and spacing constants
    panelTitleSpace = 25;
    contentHeight = 150 - panelTitleSpace;
    rowHeight = 20; % Height of each checkbox
    rowSpacing = 30; % Space between rows
    totalContentHeight = numRows * rowSpacing - (rowSpacing - rowHeight);
    
    % Calculate top position to center content vertically
    topMargin = (contentHeight - totalContentHeight) / 2;
    topCheckboxPosition = contentHeight - topMargin - rowHeight/2;
    
    % Create checkboxes in a grid layout with perfect vertical centering
    methodCheckboxes = [];
    for i = 1:numMethods
        row = ceil(i / 2);
        col = mod(i-1, 2);
        
        % Calculate positions with consistent spacing
        x = 20 + col * 340;
        y = topCheckboxPosition - (row - 1) * rowSpacing - 5;
        
        cb = uicheckbox(designPanel, 'Text', methods{i, 1}, ...
            'Position', [x y 320 20], 'Value', methods{i, 2}, ...
            'FontSize', 12);  
        
        if ~methods{i, 2}
            cb.Enable = 'off';
            cb.Value = false;
            cb.Text = [methods{i, 1}, ' (not applicable)'];
            cb.FontColor = [0.5 0.5 0.5];
        end
        
        methodCheckboxes = [methodCheckboxes; cb];
    end
    
    methodbuttony = topCheckboxPosition - 10;

    % Position Method Information button aligned with first row
    infoBtn = uibutton(designPanel, 'Text', 'Method Information', ...
        'Position', [555 methodbuttony 115 30], 'ButtonPushedFcn', @showMethodInfo, ...
        'BackgroundColor', appColors.buttonPrimary, 'FontColor', appColors.lightText);
    
    % Controller Configuration panel - perfect vertical centering
    % Calculate exact center position for perfect alignment
    panelTitleSpace = 25;
    contentHeight = 60 - panelTitleSpace;
    elementHeight = 22;
    
    % Position elements exactly in the vertical center
    verticalCenter = contentHeight/2 - elementHeight/2;
    
    structureLabel = uilabel(controllerPanel, 'Text', 'Controller Structure:', ...
        'Position', [20 verticalCenter 150 elementHeight], 'FontSize', 12);
    structureDropdown = uidropdown(controllerPanel, ...
        'Items', {'P', 'PI', 'PD', 'PID'}, ...
        'Position', [170 verticalCenter 150 elementHeight], 'Value', 'PID', ...
        'FontSize', 12);
    
    filterLabel = uilabel(controllerPanel, 'Text', 'D-Term Filter (epsilon):', ...
        'Position', [350 verticalCenter 170 elementHeight], 'FontSize', 12);
    filterField = uieditfield(controllerPanel, 'numeric', ...
        'Position', [520 verticalCenter 100 elementHeight], 'Value', 0.1, ...
        'Limits', [0.001 0.5], 'LowerLimitInclusive', true, 'UpperLimitInclusive', true, ...
        'FontSize', 12);
    
    % Performance Requirements panel - consistent spacing between elements
    panelTitleSpace = 25;
    contentHeight = 170 - panelTitleSpace;
    
    % Calculate consistent spacing
    numRows = 4;
    rowHeight = 22;
    
    % Calculate total content height and spacing to ensure equal margins
    totalRowsHeight = numRows * rowHeight;
    availableSpace = contentHeight - totalRowsHeight;
    
    % Distribute available space evenly (top margin, between rows, bottom margin)
    margin = availableSpace / (numRows + 1);
    
    % Calculate row positions with consistent spacing
    row1Y = contentHeight - margin - rowHeight;           % Top row (Damping)
    row2Y = row1Y - margin - rowHeight;                   % Second row (Phase/Robustness)
    row3Y = row2Y - margin - rowHeight;                   % Third row (Gain/Overshoot)
    row4Y = row3Y - margin - rowHeight;                   % Bottom row (Bandwidth/Goal)
    
    % Top row - Damping Ratio
    dampingLabel = uilabel(performancePanel, 'Text', 'Damping Ratio:', ...
        'Position', [20, row1Y, 150, rowHeight], 'FontSize', 12);
    dampingField = uieditfield(performancePanel, 'numeric', ...
        'Position', [170, row1Y, 100, rowHeight], 'Value', 0.8, ...
        'Limits', [0.1, 2.0], 'LowerLimitInclusive', true, 'UpperLimitInclusive', true, ...
        'FontSize', 12);
    
    % Second row - Phase Margin & Robustness
    phaseMarginLabel = uilabel(performancePanel, 'Text', 'Phase Margin (deg):', ...
        'Position', [20 row2Y 150 rowHeight], 'FontSize', 12);
    phaseMarginField = uieditfield(performancePanel, 'numeric', ...
        'Position', [170 row2Y 100 rowHeight], 'Value', 45, ...
        'Limits', [10 80], 'LowerLimitInclusive', true, 'UpperLimitInclusive', true, ...
        'FontSize', 12);
    roby = row2Y + margin + 5;
    robustnessLabel = uilabel(performancePanel, 'Text', 'Robustness:', ...
        'Position', [350 roby 100 rowHeight], 'FontSize', 12);
    robustnessDropdown = uidropdown(performancePanel, ...
        'Items', {'Low', 'Medium', 'High'}, ...
        'Position', [450 roby 100 rowHeight], 'Value', 'Medium', ...
        'FontSize', 12);
    
    % Third row - Gain Margin & Overshoot
    gainMarginLabel = uilabel(performancePanel, 'Text', 'Gain Margin (dB):', ...
        'Position', [20 row3Y 150 rowHeight], 'FontSize', 12);
    gainMarginField = uieditfield(performancePanel, 'numeric', ...
        'Position', [170 row3Y 100 rowHeight], 'Value', 8, ...
        'Limits', [3 20], 'LowerLimitInclusive', true, 'UpperLimitInclusive', true, ...
        'FontSize', 12);
    ovey = row3Y + margin + 5;
    overshootLabel = uilabel(performancePanel, 'Text', 'Overshoot (%):', ...
        'Position', [350 ovey 100 rowHeight], 'FontSize', 12);
    overshootField = uieditfield(performancePanel, 'numeric', ...
        'Position', [450 ovey 100 rowHeight], 'Value', 10, ...
        'Limits', [0 50], 'LowerLimitInclusive', true, 'UpperLimitInclusive', true, ...
        'FontSize', 12);
    
    % Fourth row - Bandwidth & Optimization Goal
    bandwidthLabel = uilabel(performancePanel, 'Text', 'Bandwidth (rad/s):', ...
        'Position', [20 row4Y 150 rowHeight], 'FontSize', 12);
    bandwidthField = uieditfield(performancePanel, 'numeric', ...
        'Position', [170 row4Y 100 rowHeight], 'Value', 1, ...
        'Limits', [0.01 100], 'LowerLimitInclusive', true, 'UpperLimitInclusive', true, ...
        'FontSize', 12);
    
    goay = row4Y + margin + 5;
    goalLabel = uilabel(performancePanel, 'Text', 'Optimization Goal:', ...
        'Position', [350 goay 120 rowHeight], 'FontSize', 12);
    goalDropdown = uidropdown(performancePanel, ...
        'Items', {'Tracking', 'Disturbance Rejection', 'Robustness'}, ...
        'Position', [450 goay 170 rowHeight], 'Value', 'Tracking', ...
        'FontSize', 12);
    
    % Score Info button aligned with the bottom row
    scoreInfoBtn = uibutton(performancePanel, 'Text', 'Score Info', ...
        'Position', [580 goay 80 rowHeight], 'ButtonPushedFcn', @showScoreInfo, ...
        'BackgroundColor', appColors.buttonPrimary, 'FontColor', appColors.lightText);
    
    % Results area with improved styling
    resultArea = uitextarea(resultsPanel, ...
        'Position', [20 60 640 90], ... % Adjusted height
        'Value', 'Auto-tuning results will be displayed here.', ...
        'Editable', 'off', 'FontSize', 12);
    
    % Buttons with improved styling and layout in a panel with no border
    buttonPanel = uipanel(resultsPanel, 'Position', [20 10 640 40], 'BorderType', 'none', 'BackgroundColor', appColors.panelBg);
    
    startBtn = uibutton(buttonPanel, 'push', 'Text', 'Start Auto-Tuning', ...
        'Position', [100 5 160 30], 'FontSize', 14, 'FontWeight', 'bold', ...
        'BackgroundColor', appColors.buttonConfirm, 'FontColor', appColors.lightText, ...
        'ButtonPushedFcn', @startAutoTuning);
    
    applyBtn = uibutton(buttonPanel, 'push', 'Text', 'Apply', ...
        'Position', [280 5 120 30], 'FontSize', 14, ...
        'Enable', 'off', 'BackgroundColor', appColors.buttonPrimary, 'FontColor', appColors.lightText, ...
        'ButtonPushedFcn', @confirmController);
    
    cancelBtn = uibutton(buttonPanel, 'push', 'Text', 'Cancel', ...
        'Position', [420 5 120 30], 'FontSize', 14, ...
        'BackgroundColor', appColors.buttonCancel, 'FontColor', appColors.lightText, ...
        'ButtonPushedFcn', @cancelDialog);
    
    % Handle figure close request
    autoTuneFig.CloseRequestFcn = @cancelDialog;
    
    % Storage for the optimized controller
    bestController = [];
    controllerDetails = '';
    
    % ========================
    % Function Implementations
    % ========================
    
    % Show method information
    function showMethodInfo(~, ~)
        methodInfoFig = uifigure('Name', 'Controller Design Methods', 'Position', [350 250 600 500]);
        methodInfoFig.Color = appColors.background;
        
        % Title panel
        methodTitlePanel = uipanel(methodInfoFig, 'Position', [10 450 580 40], 'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
        methodTitleLabel = uilabel(methodTitlePanel, 'Text', 'Controller Design Methods', ...
            'Position', [0 0 580 40], 'FontSize', 16, 'FontWeight', 'bold', ...
            'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');
        
        methodInfo = uitextarea(methodInfoFig, 'Position', [20 60 560 380], 'Editable', 'off');
        methodInfo.Value = {
            'Controller Design Methods:', 
            '-------------------------',
            '1. Ziegler-Nichols (Oscillation): Uses critical gain and frequency to determine controller parameters.',
            '   • Best for: Systems that can be safely brought to the stability boundary for testing.',
            '   • Limitations: Can result in aggressive control with high overshoot.',
            '',
            '2. Ziegler-Nichols (Step): Uses open-loop step response characteristics to determine parameters.',
            '   • Best for: Stable plants with S-shaped step response.',
            '   • Limitations: Not suitable for integrating or unstable plants.',
            '',
            '3. Aström: Improved version of Ziegler-Nichols with better robustness.',
            '   • Best for: Similar to Ziegler-Nichols, but when less overshoot is desired.',
            '',
            '4. Loop-Shaping: Frequency-domain method focused on achieving desired loop shape.',
            '   • Best for: When specific bandwidth and phase margin are required.',
            '',
            '5. IMC (Internal Model Control): Based on process model and desired closed-loop response.',
            '   • Best for: Good disturbance rejection with specified settling time.',
            '',
            '6. MIGO: Maximizes integral gain while satisfying robustness constraints.',
            '   • Best for: Good balance between performance and robustness.',
            '',
            '7. Pole Placement: Places closed-loop poles at desired locations.',
            '   • Best for: Achieving specific time domain performance.',
            '   • Features: Direct control over system dynamics.',
            '',
            '8. Compensation Controller: Directly compensates for plant dynamics.',
            '   • Best for: Systems with problematic poles or zeros.',
            '   • Features: Cancels or modifies critical plant dynamics.',
            '   • Advantages: Can handle both stable and unstable systems with proper tuning.'
        };
        
        closeBtn = uibutton(methodInfoFig, 'push', 'Text', 'Close', ...
            'Position', [250 20 100 30], 'FontSize', 12, ...
            'BackgroundColor', appColors.buttonPrimary, 'FontColor', appColors.lightText, ...
            'ButtonPushedFcn', @(~,~) close(methodInfoFig));
    end
    
    % Show score information
    function showScoreInfo(~, ~)
        scoreInfoFig = uifigure('Name', 'Controller Scoring System', 'Position', [350 300 500 400]);
        scoreInfoFig.Color = appColors.background;
        
        % Title panel
        scoreTitlePanel = uipanel(scoreInfoFig, 'Position', [10 350 480 40], 'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
        scoreTitleLabel = uilabel(scoreTitlePanel, 'Text', 'Controller Scoring System', ...
            'Position', [0 0 480 40], 'FontSize', 16, 'FontWeight', 'bold', ...
            'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');
        
        scoreInfo = uitextarea(scoreInfoFig, 'Position', [20 60 460 280], 'Editable', 'off');
        scoreInfo.Value = {
            'Controller Scoring System:', 
            '-------------------------',
            'The auto-tuning process evaluates controllers using a point system (0-100):',
            '',
            '• Excellent: 80+ points',
            '  - Meets all design requirements with excellent balance',
            '  - Optimal time and frequency domain performance',
            '  - Strong robustness to disturbances and model uncertainties',
            '',
            '• Good: 60-80 points',
            '  - Meets all key requirements with good compromise',
            '  - Stable with good time response characteristics',
            '  - Suitable phase and gain margins (>45° and >6dB)',
            '',
            '• Acceptable: 40-60 points',
            '  - Meets minimum requirements but with compromises',
            '  - May have some overshoot or longer settling time',
            '  - Sufficient but not optimal robustness',
            '',
            '• Poor: 20-40 points',
            '  - Meets basic stability requirements but performance is lacking',
            '  - May have excessive overshoot or very slow response',
            '  - Limited robustness to disturbances or model variations',
            '',
            '• Unacceptable: <20 points',
            '  - May be stable but with poor performance metrics',
            '  - Fails to meet several key design requirements',
            '  - Not recommended for implementation',
            '',
            'Note: A negative score indicates an unstable controller that should not be used.'
        };
        
        closeBtn = uibutton(scoreInfoFig, 'push', 'Text', 'Close', ...
            'Position', [200 20 100 30], 'FontSize', 12, ...
            'BackgroundColor', appColors.buttonPrimary, 'FontColor', appColors.lightText, ...
            'ButtonPushedFcn', @(~,~) close(scoreInfoFig));
    end
    
    % Function to check if a design method is applicable to the plant
    function applicable = isApplicable(G, methodName)
        % Get plant information
        [z, p, k] = zpkdata(G, 'v');
        isStable = all(real(p) < 0);
        hasIntegrator = any(abs(p) < 1e-6);
        hasRHPZero = any(real(z) > 0);
        plantOrder = length(p);
        
        % Check applicability based on plant characteristics
        switch methodName
            case 'ZN-Oscillation'
                % Works for most plants that can be stabilized
                applicable = true;
            case 'ZN-Step'
                % Only for stable plants without integrators
                applicable = isStable && ~hasIntegrator;
            case 'Astrom'
                % For stable plants
                applicable = isStable;
            case 'Loop-Shaping'
                % All plants
                applicable = true;
            case 'IMC'
                % Best for stable, minimum-phase plants but can be adapted
                applicable = true;
            case 'MIGO'
                % For stable plants primarily
                applicable = true;
            case 'Pole-Placement'
                % Pole placement is applicable to almost all systems
                applicable = true;
            case 'Compensation-Controller'
                % Compensation controller is applicable to almost all systems
                applicable = true;
            otherwise
                applicable = false;
        end
    end
    
    % Start Auto-Tuning
    function startAutoTuning(~, ~)
        % Disable the Start button during calculation
        startBtn.Enable = 'off';
        startBtn.Text = 'Processing...';
        drawnow;
        
        % Parameters structure
        options = struct();
        options.epsilon = filterField.Value;
        options.phaseMargin = phaseMarginField.Value;
        options.gainMargin = gainMarginField.Value;
        options.bandwidth = bandwidthField.Value;
        options.settlingTime = 5;  % Default value
        options.robustness = robustnessDropdown.Value;
        options.overshoot = overshootField.Value;
        options.goal = goalDropdown.Value;
        options.damping = dampingField.Value;  % Added for pole placement
        options.userSetDamping = true; % Flag to indicate user-set damping

        % Controller structure
        structure = structureDropdown.Value;
        
        % Collect selected methods
        selectedMethods = {};
        for i = 1:length(methodCheckboxes)
            if methodCheckboxes(i).Value && methodCheckboxes(i).Enable == 'on'
                methodName = methodCheckboxes(i).Text;
                % Remove "(not applicable)" if present
                methodName = strrep(methodName, ' (not applicable)', '');
                selectedMethods{end+1} = methodName;
            end
        end
        
        if isempty(selectedMethods)
            uialert(autoTuneFig, 'Please select at least one design method.', 'No Method Selected');
            startBtn.Enable = 'on';
            startBtn.Text = 'Start Auto-Tuning';
            return;
        end
        
        % Update status
        resultArea.Value = 'Auto-tuning in progress... This may take a moment.';
        drawnow;
        
        % Storage for best controller
        bestScore = -Inf;
        bestController = [];
        controllerDetails = '';
        
        try
            % Process each method
            for i = 1:length(selectedMethods)
                methodName = selectedMethods{i};
                resultArea.Value = ['Testing method: ' methodName '...'];
                drawnow;
                
                try
                    % Make sure design_controller_auto function exists and is working
                    % This is a check for the function, remove if you're sure it exists
                    if ~exist('design_controller_auto', 'file')
                        error('Function design_controller_auto not found. Please make sure it is in your MATLAB path.');
                    end
                    
                    % Map UI method name to function method name if needed
                    methodNameForFunction = methodName;
                    
                    % Design controller using selected method
                    [K_candidate, details] = design_controller_auto(G, methodNameForFunction, structure, options);
                    
                    % Extract score from details
                    score = extractScoreFromDetails(details);
                    
                    % Update the result area with current method results
                    resultArea.Value{end+1} = sprintf('Method: %s, Score: %.2f', methodName, score);
                    drawnow;
                    
                    % Check if this is the best controller so far
                    if score > bestScore
                        bestScore = score;
                        bestController = K_candidate;
                        scoreLabel = getScoreLabel(score);
                        controllerDetails = sprintf('Method: %s\nScore: %.2f (%s)\n\n%s', methodName, score, scoreLabel, details);
                    end
                catch ME
                    % Log error for this method with detailed stack trace
                    disp(['Error with method ' methodName ': ' ME.message]);
                    disp(getReport(ME, 'extended'));
                    resultArea.Value{end+1} = sprintf('Error with %s: %s', methodName, ME.message);
                    drawnow;
                end
            end
            
            % Display final results with enhanced diagnostics
            if ~isempty(bestController) && bestScore > 30
                scoreLabel = getScoreLabel(bestScore);
                resultArea.Value = {
                    'Auto-tuning completed!',
                    sprintf('Best controller designed using: %s', extractMethodName(controllerDetails)),
                    sprintf('Score: %.2f (%s)', bestScore, scoreLabel),
                    '',
                    'Check console for detailed results or click "Apply" to use this controller.'
                };
                
                % Enable Apply button
                applyBtn.Enable = 'on';
                
                % Display the controller in the command window
                disp('====== Auto-Tuning Results ======');
                disp(controllerDetails);
                
                % Get controller details for display
                [num, den] = tfdata(bestController, 'v');
                resultArea.Value{end+1} = '';
                resultArea.Value{end+1} = 'Controller Transfer Function:';
                resultArea.Value{end+1} = sprintf('Numerator: [%s]', num2str(num, '%.4g '));
                resultArea.Value{end+1} = sprintf('Denominator: [%s]', num2str(den, '%.4g '));
            elseif ~isempty(bestController) && bestScore <= 30
                % Poor result with diagnostic info
                scoreLabel = getScoreLabel(bestScore);
                
                % Get plant diagnostics
                diagnosticInfo = analyzePlantForTroubleshooting(G);
                
                resultArea.Value = {
                    'Auto-tuning completed with poor results.',
                    sprintf('Best controller score: %.2f (%s)', bestScore, scoreLabel),
                    '',
                    '== PLANT DIAGNOSTIC INFORMATION =='
                };
                
                % Add diagnostic information
                for i = 1:length(diagnosticInfo)
                    resultArea.Value{end+1} = diagnosticInfo{i};
                end
                
                resultArea.Value{end+1} = '';
                resultArea.Value{end+1} = 'RECOMMENDATIONS:';
                
                % Add recommendations based on diagnostics
                if any(real(pole(G)) > 0)
                    resultArea.Value{end+1} = '• Try the Compensation Controller for unstable systems';
                end
                
                if any(abs(pole(G)) < 1e-6)
                    resultArea.Value{end+1} = '• For systems with integrators, try PI or PID with low Ki';
                end
                
                if any(real(zero(G)) > 0)
                    resultArea.Value{end+1} = '• For non-minimum phase systems, try Compensation Controller';
                end
                
                % Recommend advanced methods
                resultArea.Value{end+1} = '• Pole Placement often works better for high-order systems';
                resultArea.Value{end+1} = '• Compensation Controller works well for multiple unstable poles';
                
                % Enable Apply button but with warning
                applyBtn.Enable = 'on';
                applyBtn.BackgroundColor = [0.9, 0.6, 0.1];  % Orange warning color
            else
                % Complete failure
                resultArea.Value = {'Auto-tuning failed. No usable controller found.'};
                
                % Get plant diagnostics
                diagnosticInfo = analyzePlantForTroubleshooting(G);
                
                resultArea.Value{end+1} = '';
                resultArea.Value{end+1} = '== REASONS FOR FAILURE ==';
                
                % Add diagnostic information
                for i = 1:length(diagnosticInfo)
                    resultArea.Value{end+1} = diagnosticInfo{i};
                end
                
                resultArea.Value{end+1} = '';
                resultArea.Value{end+1} = 'RECOMMENDATIONS:';
                resultArea.Value{end+1} = '• Try the Compensation Controller for challenging systems';
                resultArea.Value{end+1} = '• Pole Placement may work better for high-order systems';
                
                % Disable Apply button
                applyBtn.Enable = 'off';
            end
        catch ME
            % Handle general errors with detailed stack trace
            disp(['General error during auto-tuning: ' ME.message]);
            disp(getReport(ME, 'extended'));
            
            resultArea.Value = {
                'Error during auto-tuning:',
                ME.message,
                '',
                'Please try different parameters or methods.'
            };
        end
        
        % Re-enable the Start button
        startBtn.Enable = 'on';
        startBtn.Text = 'Start Auto-Tuning';
    end
    
    % Extract score from controller details
    function score = extractScoreFromDetails(details)
        % Default score if not found
        score = 0;
        
        % Look for the score information in the details string
        try
            matches = regexp(details, 'Controller Score: (\d+\.\d+)/100', 'tokens');
            if ~isempty(matches) && ~isempty(matches{1})
                score = str2double(matches{1}{1});
            else
                % If standard score format not found, look for pole placement success indicators
                if contains(details, 'Closed-loop system is stable!')
                    % For pole placement, assign a reasonable score based on stability
                    score = 75;  % Good score for a stable system
                    
                    % Reduce score if warnings are present
                    if contains(details, 'WARNING')
                        score = score - 10;
                    end
                    
                    % Boost score if it matches bandwidth and damping requirements
                    if contains(details, 'Successfully computed state feedback')
                        score = score + 5;
                    end
                else
                    % If system is unstable, assign a low score
                    score = 20;  % Poor score
                end
            end
        catch
            % If there's any error, return default score
            score = 0;
        end
    end
    
    % Helper function to extract method name from details
    function methodName = extractMethodName(details)
        lines = strsplit(details, '\n');
        if ~isempty(lines) && length(lines) > 0
            firstLine = lines{1};
            parts = strsplit(firstLine, ':');
            if length(parts) > 1
                methodName = strtrim(parts{2});
            else
                methodName = 'Unknown';
            end
        else
            methodName = 'Unknown';
        end
    end
    
    % Helper function to get a textual label for a score
    function label = getScoreLabel(score)
        if score < 0
            label = 'Unstable';
        elseif score < 20
            label = 'Unacceptable';
        elseif score < 40
            label = 'Poor';
        elseif score < 60
            label = 'Acceptable';
        elseif score < 80
            label = 'Good';
        else
            label = 'Excellent';
        end
    end
    
    % Confirm Controller Function
    function confirmController(~, ~)
        if ~isempty(bestController)
            % Set K to the best controller and close dialog
            K = bestController;
            uiresume(mainFig);
            close(autoTuneFig);
            close(mainFig);
        else
            uialert(autoTuneFig, 'No valid controller available. Please run auto-tuning first.', 'No Controller');
        end
    end
    
    % Cancel Function
    function cancelDialog(~, ~)
        % Close auto-tuning window without selecting a controller
        close(autoTuneFig);
    end
    
    % Plant diagnostics function
    function diagnosticInfo = analyzePlantForTroubleshooting(G)
        % Analyzes a plant transfer function for troubleshooting
        % Returns a cell array of diagnostic messages
        
        diagnosticInfo = {};
        
        try
            % Get poles and zeros
            p = pole(G);
            try
                z = zero(G);
            catch
                z = [];
            end
            
            % Get transfer function data
            [num, den] = tfdata(G, 'v');
            
            % Check system order
            order = length(p);
            if order > 3
                diagnosticInfo{end+1} = sprintf('• High-order system (order %d)', order);
            end
            
            % Check stability
            unstable_poles = p(real(p) > 0);
            if ~isempty(unstable_poles)
                diagnosticInfo{end+1} = sprintf('• Unstable system with %d RHP pole(s)', length(unstable_poles));
                for i = 1:min(length(unstable_poles), 3)
                    if imag(unstable_poles(i)) ~= 0
                        diagnosticInfo{end+1} = sprintf('  - Complex pole at %.3f + %.3fi', real(unstable_poles(i)), imag(unstable_poles(i)));
                    else
                        diagnosticInfo{end+1} = sprintf('  - Real pole at %.3f', real(unstable_poles(i)));
                    end
                end
            end
            
            % Check for integrators
            integrators = p(abs(p) < 1e-6);
            if ~isempty(integrators)
                diagnosticInfo{end+1} = sprintf('• System has %d integrator(s)', length(integrators));
            end
            
            % Check for RHP zeros
            if ~isempty(z)
                rhp_zeros = z(real(z) > 0);
                if ~isempty(rhp_zeros)
                    diagnosticInfo{end+1} = sprintf('• Non-minimum phase with %d RHP zero(s)', length(rhp_zeros));
                end
            end
            
            % Check for proper/strictly proper
            rel_degree = length(den) - length(num);
            if rel_degree < 0
                diagnosticInfo{end+1} = '• Improper transfer function (numerator order > denominator order)';
            elseif rel_degree == 0
                diagnosticInfo{end+1} = '• Proper but not strictly proper transfer function';
            end
            
            % Check for oscillatory modes
            osc_poles = p(abs(imag(p)) > 0.1*abs(real(p)) & real(p) < 0);
            if ~isempty(osc_poles)
                diagnosticInfo{end+1} = sprintf('• System has %d oscillatory mode(s)', length(osc_poles)/2);
            end
            
            % Check for poorly damped modes
            poor_damp = false;
            for i = 1:length(p)
                if real(p(i)) < 0 && imag(p(i)) ~= 0
                    damp_ratio = -real(p(i))/abs(p(i));
                    if damp_ratio < 0.2
                        poor_damp = true;
                        break;
                    end
                end
            end
            if poor_damp
                diagnosticInfo{end+1} = '• System has poorly damped modes (damping ratio < 0.2)';
            end
            
            % Special cases for auto-tuning failure
            if ~isempty(unstable_poles) && length(unstable_poles) > 1
                diagnosticInfo{end+1} = '• Multiple unstable poles make classical tuning methods fail';
            end
            
            if ~isempty(rhp_zeros) && ~isempty(unstable_poles)
                diagnosticInfo{end+1} = '• Combination of RHP zeros and poles creates fundamental limitations';
            end
            
            if rel_degree <= 0 && ~isempty(unstable_poles)
                diagnosticInfo{end+1} = '• Improper/non-strictly proper transfer function with instability';
                diagnosticInfo{end+1} = '  This combination is particularly challenging for control design';
            end
            
            % Add overall assessment
            if ~isempty(unstable_poles) || rel_degree <= 0 || (~isempty(z) && ~isempty(rhp_zeros))
                diagnosticInfo{end+1} = '';
                diagnosticInfo{end+1} = 'OVERALL: This system has challenging characteristics';
                diagnosticInfo{end+1} = 'Consider using the Compensation Controller design method';
            end
            
        catch ME
            % If analysis fails, return basic error info
            diagnosticInfo{end+1} = 'Error during plant analysis:';
            diagnosticInfo{end+1} = ['• ' ME.message];
        end
    end
end

    % Helper function: Computes the controller transfer function based on controller type and inputs
    function Ktemp = computeController(ctrlType, fields)
        switch ctrlType
            case 'P'
                Kp = convertValue(fields{1}.Value);
                Ktemp = tf(Kp, 1);
            case 'PI'
                Kp = convertValue(fields{1}.Value);
                Ki = convertValue(fields{2}.Value);
                Ktemp = tf([Kp, Ki], [1, 0]);
            case 'PD'
                Kp = convertValue(fields{1}.Value);
                Kd = convertValue(fields{2}.Value);
                % Apply filter for derivative term (0.1 = default epsilon)
                Ktemp = tf([Kd, Kp], [0.1*Kd, 1]);
            case 'PID'
                Kp = convertValue(fields{1}.Value);
                Ki = convertValue(fields{2}.Value);
                Kd = convertValue(fields{3}.Value);
                % Apply filter for derivative term (0.1 = default epsilon)
                Ktemp = tf([Kd, Kp, Ki], [0.1*Kd, 1, 0]);
            case 'PT1'
                Kp = convertValue(fields{1}.Value);
                T  = convertValue(fields{2}.Value);
                Ktemp = tf(Kp, [T, 1]);
            case 'PIT1'
                Kp = convertValue(fields{1}.Value);
                Ki = convertValue(fields{2}.Value);
                T  = convertValue(fields{3}.Value);
                Ktemp = tf([Kp, Ki], [T, 1, 0]);
            case 'I2'
                Ki = convertValue(fields{1}.Value);
                Ktemp = tf(Ki, [1, 0, 0]);
            case 'PIDT1'
                Kp = convertValue(fields{1}.Value);
                Ki = convertValue(fields{2}.Value);
                Kd = convertValue(fields{3}.Value);
                T  = convertValue(fields{4}.Value);
                Ktemp = tf([Kd, Kp, Ki], [T, 1, 0]);
            case 'Custom'
                numStr = strrep(fields{1}.Value, ',', '.');
                denStr = strrep(fields{2}.Value, ',', '.');
                numVec = str2num(numStr);  %#ok<ST2NM>
                denVec = str2num(denStr);  %#ok<ST2NM>
                if isempty(numVec) || isempty(denVec) || any(isnan(numVec)) || any(isnan(denVec))
                    error('Invalid input');
                end
                Ktemp = tf(numVec, denVec);
            case 'Pole Placement'
                % Extract parameters
                bandwidth = convertValue(fields{1}.Value);
                damping = convertValue(fields{2}.Value);
                epsilon = convertValue(fields{3}.Value);
                
                % Set up options for designPolePlacement
                options = struct('bandwidth', bandwidth, 'damping', damping, 'epsilon', epsilon);
                
                % Perform plant analysis
                plantInfo = analyzePlant(G);
                
                try
                    % Try to use designPolePlacement function if available
                    if exist('designPolePlacement', 'file')
                        [Ktemp, ~] = designPolePlacement(G, 'PID', options, plantInfo);
                    else
                        % If function not found, create a basic PID directly
                        error('designPolePlacement function not found');
                    end
                catch ME
                    % Log the error
                    disp(['Pole placement error: ' ME.message]);
                    disp('Creating fallback controller');
                    
                    % Create a conservative PID controller directly
                    if plantInfo.isUnstable
                        % For unstable plants with poles in the right half plane
                        p = plantInfo.poles;
                        unstable_poles = p(real(p) > 0);
                        
                        if ~isempty(unstable_poles)
                            % Find the rightmost pole (most unstable)
                            [~, idx] = max(real(unstable_poles));
                            max_real_part = real(unstable_poles(idx));
                            
                            % Create a stabilizing controller
                            Kp = max_real_part * 2;
                            Ki = 0.01;
                            Kd = Kp * 5;
                        else
                            % Default values if we can't determine unstable poles
                            Kp = 1;
                            Ki = 0.01;
                            Kd = 5;
                        end
                    else
                        % For stable plants, use conservative values
                        Kp = 0.5;
                        Ki = 0.05;
                        Kd = 2;
                    end
                    
                    % Create PID controller with derivative filtering
                    Ktemp = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                    
                    % Issue a warning in the command window
                    warning('Pole placement failed. Using fallback PID controller.');
                end
            otherwise
                error('Unknown controller type.');
        end
    end

    % Helper function: Converts an input string (allowing decimal commas) to a number.
    % If the input string contains spaces, it takes only the first token.
    function num = convertValue(str)
        % Replace commas with periods
        str = strrep(str, ',', '.');
        % Tokenize the string at whitespace and take only the first token
        token = strtok(str);
        num = str2double(token);
        if isnan(num)
            error('Conversion error');
        end
    end

    % Helper function: Converts a polynomial coefficient vector to an HTML-formatted
    % polynomial string without an asterisk between the coefficient and 's', using <sup> tags for exponents.
    function polyStr = polyToHTMLString(coeff)
        deg = length(coeff) - 1;  % Determine the degree of the polynomial
        terms = {};                % Cell array to store individual terms

        % Iterate over all coefficients
        for i = 1:length(coeff)
            coef = coeff(i);
            exp = deg - (i - 1);   % Determine the exponent for the current term

            % Skip zero coefficients as they don't affect the expression
            if coef == 0
                continue;
            end

            % Handle signs: add " + " if not the first term;
            % For negative coefficients, add " - " and use the absolute value.
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

            % Add "s" and superscript for exponents greater than 0.
            if exp > 1
                term = strcat(term, 's<sup>', num2str(exp), '</sup>');
            elseif exp == 1
                term = strcat(term, 's');
            end

            % Add the formatted term to the cell array.
            terms{end + 1} = term;
        end

        % Join all terms into a single string
        polyStr = strjoin(terms, '');
        
        % If all coefficients are zero, return "0"
        if isempty(polyStr)
            polyStr = '0';
        end
    end

    % Helper function: Returns an HTML string with the symbolic preview
    function htmlStr = getSymbolicPreviewHTML(ctrlType, ~)
        
        switch ctrlType
            case 'P'
                numHTML = 'Kp';
                denHTML = '1';
            case 'PI'
                numHTML = 'Kp s + Ki';
                denHTML = 's';
            case 'PD'
                numHTML = 'Kp + Kd s';
                denHTML = '1';
            case 'PID'
                numHTML = 'Kd s<sup>2</sup> + Kp s + Ki';
                denHTML = 's';
            case 'PT1'
                numHTML = 'Kp';
                denHTML = 'T s + 1';
            case 'PIT1'
                numHTML = 'Kp Ki s + Kp';
                denHTML = 'T s<sup>2</sup> + s';
            case 'I2'
                numHTML = 'Ki';
                denHTML = 's<sup>2</sup>';
            case 'PIDT1'
                numHTML = 'Kd s<sup>2</sup> + Kp s + Ki';
                denHTML = 'T s<sup>2</sup> + s';
            case 'Custom'
                numHTML = 'numerator';
                denHTML = 'denominator';
            case 'Pole Placement'
                numHTML = 'Controller based on desired poles';
                denHTML = 'with bandwidth ω and damping ζ';
            otherwise
                numHTML = '?';
                denHTML = '?';
        end
        
        % Create the HTML fraction with exact same styling as transfer function UI
        htmlStr = ['<html><head><style>', ...
            'body { display: flex; justify-content: center; align-items: center; height: 100%; margin: 0; padding: 0; }', ...
            '.fraction { display: inline-block; vertical-align: middle; margin: 0 auto; text-align: center; }', ...
            '.fraction .num { border-bottom: 1px solid black; padding: 5px 15px; font-size: 16px; }', ...
            '.fraction .den { padding: 5px 15px; font-size: 16px; }', ...
            '</style></head><body>', ...
            '<div class="fraction">', ...
            '<div class="num">', numHTML, '</div>', ...
            '<div class="den">', denHTML, '</div>', ...
            '</div>', ...
            '</body></html>'];
    end

    % Plant analysis function for the pole placement controller
    function plantInfo = analyzePlant(G)
        % ANALYZEPLANT Analyze plant characteristics to guide controller design
        % Outputs a structure with plant information
        
        plantInfo = struct();
        
        % Get poles and zeros
        plantInfo.poles = pole(G);
        try
            plantInfo.zeros = zero(G);
        catch
            plantInfo.zeros = [];
        end
        
        % Check stability
        plantInfo.isUnstable = any(real(plantInfo.poles) > 0);
        
        % Check for integrators (poles at the origin)
        plantInfo.hasIntegrator = any(abs(plantInfo.poles) < 1e-6);
        
        % Check for non-minimum phase zeros (RHP zeros)
        plantInfo.hasRHPZeros = any(real(plantInfo.zeros) > 0);
        
        % Check for time delay
        [num, den] = tfdata(G, 'v');
        % Pade approximations typically have alternating sign coefficients
        if length(num) > 1 && all(sign(num(1:2:end)) ~= sign(num(2:2:end)))
            plantInfo.hasDelay = true;
        else
            plantInfo.hasDelay = false;
        end
        
        % Determine plant DC gain
        try
            plantInfo.dcGain = dcgain(G);
        catch
            % For plants with pure integrators
            plantInfo.dcGain = Inf;
        end
        
        % Check if high order (more than 2 states)
        plantInfo.isHighOrder = length(plantInfo.poles) > 2;
        
        % Get step response characteristics if plant is stable
        if ~plantInfo.isUnstable
            try
                t = linspace(0, 100, 1000);
                [y, t] = step(G, t);
                plantInfo.stepInfo = stepinfo(y, t);
                plantInfo.stepResponse = struct('time', t, 'response', y);
                
                % Add FOPDT model fields with empty values
                plantInfo.FOPDT = struct('L', NaN, 'T', NaN, 'K', NaN);
            catch
                plantInfo.stepInfo = [];
                plantInfo.stepResponse = [];
                plantInfo.FOPDT = struct('L', NaN, 'T', NaN, 'K', NaN);
            end
        else
            plantInfo.stepInfo = [];
            plantInfo.stepResponse = [];
            plantInfo.FOPDT = struct('L', NaN, 'T', NaN, 'K', NaN);
        end
        
        % For unstable plants, estimate stabilizing feedback
        if plantInfo.isUnstable
            % Simple estimate for stabilizing gain
            p = plantInfo.poles;
            [maxRealPart, idx] = max(real(p));
            if maxRealPart <= 0
                plantInfo.stabilizingGain = 0;
            elseif imag(p(idx)) == 0
                plantInfo.stabilizingGain = maxRealPart * 1.5;
            else
                plantInfo.stabilizingGain = maxRealPart * 2;
            end
        else
            plantInfo.stabilizingGain = 0;
        end
    end
end