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

% Create UI figure with enhanced styling
cont_fig = uifigure('Name', 'Reference Signal Input', 'Position', [400, 200, 600, 500]);
cont_fig.Color = appColors.background;

% Add title panel with enhanced styling
titlePanel = uipanel(cont_fig, 'Position', [10 450 580 40], 'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
titleLabel = uilabel(titlePanel, 'Text', 'Continuous Function Input (r(t))', ...
    'Position', [0 0 580 40], 'FontSize', 16, 'FontWeight', 'bold', ...
    'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');

% Add plot area in a panel
previewPanel = uipanel(cont_fig, 'Title', 'Function Preview', ...
    'Position', [10 250 580 190], 'TitlePosition', 'centertop', ...
    'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

ax = uiaxes(previewPanel, 'Position', [20, 20, 540, 140]);
xlabel(ax, 'Time (s)');
ylabel(ax, 'r(t)');
title(ax, 'Preview');
grid(ax, 'on');

% Time parameters panel
timePanel = uipanel(cont_fig, 'Title', 'Time Parameters', ...
    'Position', [10 170 580 70], 'TitlePosition', 'centertop', ...
    'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

% Time vector inputs with improved layout
uilabel(timePanel, 'Text', 'Start Time:', 'Position', [20, 20, 80, 22], 'FontSize', 12);
startField = uieditfield(timePanel, 'numeric', 'Position', [110, 20, 80, 22], 'Value', 0, 'FontSize', 12);

uilabel(timePanel, 'Text', 'End Time:', 'Position', [210, 20, 80, 22], 'FontSize', 12);
endField = uieditfield(timePanel, 'numeric', 'Position', [290, 20, 80, 22], 'Value', 10, 'FontSize', 12);

uilabel(timePanel, 'Text', 'Time Steps:', 'Position', [390, 20, 80, 22], 'FontSize', 12);
stepField = uieditfield(timePanel, 'numeric', 'Position', [480, 20, 80, 22], 'Value', 0.1, 'FontSize', 12);

% Function input panel
functionPanel = uipanel(cont_fig, 'Title', 'Function Definition', ...
    'Position', [10 100 580 60], 'TitlePosition', 'centertop', ...
    'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

uilabel(functionPanel, 'Text', 'r(t) =', 'Position', [20, 15, 50, 22], 'FontSize', 12, 'FontWeight', 'bold');
funcField = uieditfield(functionPanel, 'text', 'Position', [70, 15, 410, 22], 'Value', '', 'FontSize', 12);

% Add help button in the function panel
helpBtn = uibutton(functionPanel, 'push', 'Text', '?', ...
    'Position', [490, 15, 30, 22], ...
    'BackgroundColor', appColors.buttonPrimary, ...
    'FontColor', appColors.lightText, ...
    'FontWeight', 'bold', ...
    'FontSize', 12, ...
    'ButtonPushedFcn', @(btn, event) showHelpDialog());

% Button panel
buttonPanel = uipanel(cont_fig, 'Position', [10 10 580 80], 'BorderType', 'none', 'BackgroundColor', appColors.background);

% Preview button
previewButton = uibutton(buttonPanel, 'push', 'Position', [140, 25, 100, 30], ...
    'Text', 'Preview', ...
    'FontSize', 12, ...
    'BackgroundColor', appColors.buttonPrimary, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) target_preview_continuous_function(ax, startField, endField, stepField, funcField));

% Confirm button
confirmButton = uibutton(buttonPanel, 'push', 'Position', [250, 25, 100, 30], ...
    'Text', 'Confirm', ...
    'FontSize', 12, ...
    'BackgroundColor', appColors.buttonConfirm, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) target_confirm_continuous_function(selection_fig, cont_fig, startField, endField, stepField, funcField));

% Cancel button
cancelButton = uibutton(buttonPanel, 'push', 'Position', [360, 25, 100, 30], ...
    'Text', 'Cancel', ...
    'FontSize', 12, ...
    'BackgroundColor', appColors.buttonCancel, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) target_cancel_input(cont_fig));

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