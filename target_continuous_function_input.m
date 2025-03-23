function target_continuous_function_input(selection_fig)
% TARGET_CONTINUOUS_FUNCTION_INPUT - Handle continuous function input for reference signal
%
% Parameters:
%   selection_fig - The main selection figure handle

% Create UI figure for continuous function input
cont_fig = uifigure('Name', 'Continuous Function Input', 'Position', [400, 300, 600, 450]);

% Add plot area at top, centered
ax = uiaxes(cont_fig, 'Position', [100, 200, 400, 200]);
xlabel(ax, 'Time (s)');
ylabel(ax, 'r(t)');
title(ax, 'Function Preview');
grid(ax, 'on');

% Time vector inputs - positioned in a row under the plot
pnl = uipanel(cont_fig, 'Position', [100, 140, 400, 30], 'BorderType', 'none');

% Start Time
uilabel(pnl, 'Text', 'Start Time:', 'Position', [10, 5, 60, 22]);
startField = uieditfield(pnl, 'numeric', 'Position', [75, 5, 50, 22], 'Value', 0);

% End Time
uilabel(pnl, 'Text', 'End Time:', 'Position', [135, 5, 60, 22]);
endField = uieditfield(pnl, 'numeric', 'Position', [195, 5, 50, 22], 'Value', 10);

% Time Steps
uilabel(pnl, 'Text', 'Time Steps:', 'Position', [255, 5, 60, 22]);
stepField = uieditfield(pnl, 'numeric', 'Position', [320, 5, 50, 22], 'Value', 0.1);

% Function input with clear label
uilabel(cont_fig, 'Text', 'r(t):', 'Position', [150, 100, 40, 22]);
funcField = uieditfield(cont_fig, 'text', 'Position', [190, 100, 260, 22], 'Value', '');

% Function info text
uilabel(cont_fig, 'Text', 'Enter any valid MATLAB expression using "t" as the variable.', ...
    'Position', [100, 70, 400, 22]);
uilabel(cont_fig, 'Text', 'Examples: sin(t), t^2, sqrt(t), pi*cos(t), exp(-t)', ...
    'Position', [100, 50, 400, 22]);

% Buttons at bottom
buttonWidth = 100;
buttonHeight = 25;
panelWidth = cont_fig.Position(3);
buttonY = 15;

% Preview button
previewButton = uibutton(cont_fig, 'Position', [(panelWidth/2)-150, buttonY, buttonWidth, buttonHeight], ...
    'Text', 'Preview', ...
    'ButtonPushedFcn', @(btn, event) target_preview_continuous_function(ax, startField, endField, stepField, funcField));

% Confirm button
confirmButton = uibutton(cont_fig, 'Position', [(panelWidth/2)-50, buttonY, buttonWidth, buttonHeight], ...
    'Text', 'Confirm', ...
    'ButtonPushedFcn', @(btn, event) target_confirm_continuous_function(selection_fig, cont_fig, startField, endField, stepField, funcField));

% Cancel button
cancelButton = uibutton(cont_fig, 'Position', [(panelWidth/2)+50, buttonY, buttonWidth, buttonHeight], ...
    'Text', 'Cancel', ...
    'ButtonPushedFcn', @(btn, event) target_cancel_input(cont_fig));

% Handle window close event
cont_fig.CloseRequestFcn = @(src, event) target_cancel_input(src);

% Wait for user to confirm or cancel
uiwait(cont_fig);
end