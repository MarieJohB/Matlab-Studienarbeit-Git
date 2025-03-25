function disturbance_continuous_function_input(selection_fig, disturbance_label)
% DISTURBANCE_CONTINUOUS_FUNCTION_INPUT - Enhanced UI for continuous disturbance function input
%
% Parameters:
%   selection_fig - The main selection figure handle
%   disturbance_label - Label for the disturbance (e.g. 'd_1')

% Define unified color scheme
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
cont_fig = uifigure('Name', [disturbance_label ' Function Input'], 'Position', [400, 200, 600, 550]);
cont_fig.Color = appColors.background;

% Add title panel with enhanced styling
titlePanel = uipanel(cont_fig, 'Position', [10 500 580 40], 'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
titleLabel = uilabel(titlePanel, 'Text', ['Continuous Function Input (', disturbance_label, '(t))'], ...
    'Position', [0 0 580 40], 'FontSize', 16, 'FontWeight', 'bold', ...
    'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');

% Add plot area in a panel - MAKE BIGGER
previewPanel = uipanel(cont_fig, 'Title', 'Function Preview', ...
    'Position', [10 240 580 250], 'TitlePosition', 'centertop', ...
    'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

ax = uiaxes(previewPanel, 'Position', [20, 20, 540, 200]);
xlabel(ax, 'Time (s)');
ylabel(ax, disturbance_label);
% No title
grid(ax, 'on');

% Time parameters panel
timePanel = uipanel(cont_fig, 'Title', 'Time Parameters', ...
    'Position', [10 160 580 70], 'TitlePosition', 'centertop', ...
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
    'Position', [10 90 580 60], 'TitlePosition', 'centertop', ...
    'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

uilabel(functionPanel, 'Text', [disturbance_label, '(t) ='], 'Position', [20, 15, 50, 22], 'FontSize', 12, 'FontWeight', 'bold');
funcField = uieditfield(functionPanel, 'text', 'Position', [70, 15, 490, 22], 'Value', '', 'FontSize', 12);

% Button panel with title
buttonPanel = uipanel(cont_fig, 'Title', 'Actions', ...
    'Position', [10 10 580 70], 'TitlePosition', 'centertop', ...
    'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

% Center buttons in panel
panelCenter = 580/2;
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
    'ButtonPushedFcn', @(btn, event) disturbance_preview_continuous_function(ax, startField, endField, stepField, funcField, disturbance_label));

% Confirm button
confirmButton = uibutton(buttonPanel, 'push', 'Position', [startX + buttonWidth + spacing, 25, buttonWidth, 30], ...
    'Text', 'Confirm', ...
    'FontSize', 12, ...
    'BackgroundColor', appColors.buttonConfirm, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) disturbance_confirm_continuous_function(selection_fig, cont_fig, startField, endField, stepField, funcField));

% Cancel button
cancelButton = uibutton(buttonPanel, 'push', 'Position', [startX + 2*(buttonWidth + spacing), 25, buttonWidth, 30], ...
    'Text', 'Cancel', ...
    'FontSize', 12, ...
    'BackgroundColor', appColors.buttonCancel, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) disturbance_cancel_input(cont_fig));

% Help button renamed and centered under Confirm button
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
        ['Function Input Guide for ' disturbance_label '(t):'], 
        '-----------------------------------------', 
        'Enter any valid MATLAB expression using "t" as the time variable.', 
        '', 
        'Examples:', 
        '• Constants: "5" or "3.14"', 
        '• Linear functions: "2*t" or "t+5"', 
        '• Polynomials: "t^2 + 3*t + 1"', 
        '• Trigonometric: "sin(t)" or "cos(2*t)"', 
        '• Exponentials: "exp(-t)" or "exp(-0.5*t)*sin(t)"', 
        '• Step functions: "t>=2" (returns 1 when t≥2, 0 otherwise)', 
        '• Random noise: "randn(1)" (standard normal distribution)', 
        '• Impulse: "(t>1 && t<1.1)"', 
        '', 
        'For disturbances, consider:', 
        '• Step inputs: "2*(t>5)" (step at t=5 with amplitude 2)', 
        '• Random disturbances: "0.5*randn(1)"', 
        '• Sinusoidal disturbances: "0.2*sin(2*t)"', 
        '• Bounded noise: "0.1*sign(randn(1))"', 
        '', 
        'Tips:', 
        '• Use "*" for multiplication: "5*t" not "5t"', 
        '• Ensure proper syntax for functions and operators', 
        '• For no disturbance, enter "0" or leave empty'
    };
    
    % Close button
    closeBtn = uibutton(helpFig, 'push', 'Text', 'Close', ...
        'Position', [200, 20, 100, 30], ...
        'BackgroundColor', appColors.buttonPrimary, ...
        'FontColor', appColors.lightText, ...
        'ButtonPushedFcn', @(~,~) close(helpFig));
end

% Handle window close event
cont_fig.CloseRequestFcn = @(src, event) disturbance_cancel_input(src);

% Wait for user to confirm or cancel
uiwait(cont_fig);
end