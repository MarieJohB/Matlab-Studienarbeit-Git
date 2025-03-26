function target_continuous_function_input(selection_fig)
% TARGET_CONTINUOUS_FUNCTION_INPUT - Enhanced UI for continuous function input
%
% Parameters:
%   selection_fig - The main selection figure handle

% Define unified color scheme to match get_user_controller and get_user_transfer_function
appColors = struct(...
    'background', [0.95 0.95 0.97], ...    % Light gray background
    'panelHeader', [0.2 0.4 0.7], ...      % Blue panel header
    'panelBg', [0.95 0.95 0.97], ...       % Light panel background
    'buttonPrimary', [0.3 0.6 0.9], ...    % Blue buttons
    'buttonConfirm', [0.3 0.8 0.3], ...    % Green confirm button
    'buttonCancel', [0.8 0.3 0.3], ...     % Red cancel button
    'buttonHelp', [0.5 0.5 0.5], ...       % Gray help button
    'text', [0.2 0.2 0.2], ...             % Dark text
    'lightText', [1 1 1]);                 % White text for dark backgrounds

% Create UI figure with enhanced styling - INCREASED SIZE
cont_fig = uifigure('Name', 'Reference Signal Input', 'Position', [400, 200, 700, 615]);
cont_fig.Color = appColors.background;

% Add title panel with enhanced styling
titlePanel = uipanel(cont_fig, 'Position', [10 565 680 40], 'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
titleLabel = uilabel(titlePanel, 'Text', 'Continuous Function Input (r(t))', ...
    'Position', [0 0 695 40], 'FontSize', 16, 'FontWeight', 'bold', ...
    'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');

% Add plot area in a panel - INCREASED SIZE
previewPanel = uipanel(cont_fig, 'Title', 'Function Preview', ...
    'Position', [10 255 680 300], 'TitlePosition', 'centertop', ...
    'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

ax = uiaxes(previewPanel, 'Position', [20, 10, 640, 250]);
xlabel(ax, 'Time (s)');
% No y-label as requested
grid(ax, 'on');

% Time parameters panel
timePanel = uipanel(cont_fig, 'Title', 'Time Parameters', ...
    'Position', [10 175 680 70], 'TitlePosition', 'centertop', ...
    'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

% Calculate center positions for time parameter fields
panelWidth = 680;
totalFieldsWidth = 680;  % Total width of all labels and fields
startX = (panelWidth - totalFieldsWidth) / 2;

% Time vector inputs with centered layout
uilabel(timePanel, 'Text', 'Start Time:', 'Position', [10, 13, 80, 22], 'FontSize', 12);
startField = uieditfield(timePanel, 'numeric', 'Position', [90, 13, 80, 22], 'Value', 0, 'FontSize', 12);

uilabel(timePanel, 'Text', 'End Time:', 'Position', [220, 13, 80, 22], 'FontSize', 12);
endField = uieditfield(timePanel, 'numeric', 'Position', [300, 13, 80, 22], 'Value', 10, 'FontSize', 12);

uilabel(timePanel, 'Text', 'Time Steps:', 'Position', [430, 13, 80, 22], 'FontSize', 12);
stepField = uieditfield(timePanel, 'numeric', 'Position', [510, 13, 80, 22], 'Value', 0.1, 'FontSize', 12);

% Function input panel
functionPanel = uipanel(cont_fig, 'Title', 'Function Definition', ...
    'Position', [10 105 680 60], 'TitlePosition', 'centertop', ...
    'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

uilabel(functionPanel, 'Text', 'r(t) =', 'Position', [20, 8, 50, 22], 'FontSize', 12, 'FontWeight', 'bold');
funcField = uieditfield(functionPanel, 'text', 'Position', [70, 8, 590, 22], 'Value', '', 'FontSize', 12);

% Button panel - INCREASED SIZE
buttonPanel = uipanel(cont_fig, 'Title', 'Actions', ...
    'Position', [10 10 680 85], 'TitlePosition', 'centertop', ...
    'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

% Center buttons in panel
panelCenter = 680/2;
buttonWidth = 100;
spacing = 20;
totalWidth = 3*buttonWidth + 2*spacing;
startX = panelCenter - totalWidth/2;

% Preview button
previewButton = uibutton(buttonPanel, 'push', 'Position', [startX, 25, buttonWidth, 30], ...
    'Text', 'Preview', ...
    'FontSize', 12, ...
    'BackgroundColor', appColors.buttonPrimary, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) target_preview_continuous_function(ax, startField, endField, stepField, funcField));

% Confirm button
confirmButton = uibutton(buttonPanel, 'push', 'Position', [startX + buttonWidth + spacing, 25, buttonWidth, 30], ...
    'Text', 'Confirm', ...
    'FontSize', 12, ...
    'BackgroundColor', appColors.buttonConfirm, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) target_confirm_continuous_function(selection_fig, cont_fig, startField, endField, stepField, funcField));

% Cancel button
cancelButton = uibutton(buttonPanel, 'push', 'Position', [startX + 2*(buttonWidth + spacing), 25, buttonWidth, 30], ...
    'Text', 'Cancel', ...
    'FontSize', 12, ...
    'BackgroundColor', appColors.buttonCancel, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) target_cancel_input(cont_fig));

% Help button centered under Confirm button
helpBtn = uibutton(buttonPanel, 'push', 'Text', 'Help', ...
    'Position', [startX + buttonWidth + spacing, 5, buttonWidth, 18], ...
    'BackgroundColor', appColors.buttonPrimary, ...
    'FontColor', appColors.lightText, ...
    'FontSize', 10, ...
    'ButtonPushedFcn', @(btn, event) showHelpDialog());

% Function to show help dialog
function showHelpDialog()
    helpFig = uifigure('Name', 'Function Input Help', 'Position', [450, 300, 500, 400]);
    helpFig.Color = appColors.background;
    
    % Title panel for help
    helpTitlePanel = uipanel(helpFig, 'Position', [10 350 480 40], ...
        'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
    helpTitleLabel = uilabel(helpTitlePanel, 'Text', 'Function Input Reference', ...
        'Position', [0 0 480 40], 'FontSize', 16, 'FontWeight', 'bold', ...
        'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');
    
    % Help text area
    helpText = uitextarea(helpFig, 'Position', [20 60 460 280], 'Editable', 'off', 'FontSize', 12);
    helpText.Value = {
        'Function Input Guide:', 
        '-----------------------', 
        'Enter any valid MATLAB expression using "t" as the time variable.', 
        '', 
        'Examples:', 
        '• Constants: "5" or "3.14"', 
        '• Linear functions: "2*t" or "t+5"', 
        '• Polynomials: "t^2 + 3*t + 1"', 
        '• Trigonometric: "sin(t)" or "cos(2*t)"', 
        '• Exponentials: "exp(-t)" or "exp(-0.5*t)*sin(t)"', 
        '• Step functions: "t>=2" (returns 1 when t≥2, 0 otherwise)', 
        '• Conditional: "(t<5)*sin(t) + (t>=5)*cos(t)"', 
        '• Damped oscillation: "exp(-0.2*t)*sin(2*t)"', 
        '', 
        'Tips:', 
        '• Use "*" for multiplication: "5*t" not "5t"', 
        '• Ensure proper syntax for functions and operators', 
        '• You can use mathematical constants like pi'
    };
    
    % Close button
    closeBtn = uibutton(helpFig, 'push', 'Text', 'Close', ...
        'Position', [200, 20, 100, 30], ...
        'BackgroundColor', appColors.buttonPrimary, ...
        'FontColor', appColors.lightText, ...
        'ButtonPushedFcn', @(~,~) close(helpFig));
end

% Handle window close event
cont_fig.CloseRequestFcn = @(src, event) target_cancel_input(src);

% Wait for user to confirm or cancel
uiwait(cont_fig);
end