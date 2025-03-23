function disturbance_piecewise_function_input(selection_fig, disturbance_label)
% DISTURBANCE_PIECEWISE_FUNCTION_INPUT - Handle piecewise function input for disturbance signal
%
% Parameters:
%   selection_fig - The main selection figure handle
%   disturbance_label - The label for the disturbance (e.g. 'd_1')

% Initialize variable for number of sections
num_sections = [];

while isempty(num_sections)
    prompt = {'Enter the number of sections (max 5): '};
    dlgtitle = 'Number of Sections';
    dims = [1 50];
    answer = inputdlg(prompt, dlgtitle, dims);

    % Check if user cancelled the dialog
    if isempty(answer)
        disp('Operation cancelled by user.');
        return;
    end

    % Replace commas with periods and parse number
    num_sections_str = strrep(answer{1}, ',', '.');
    num_sections = str2double(num_sections_str);

    % Validate input
    if isnan(num_sections) || num_sections < 1 || num_sections > 5 || mod(num_sections, 1) ~= 0
        uiwait(msgbox('Invalid input. Please enter a number between 1 and 5.', 'Error', 'error'));
        num_sections = [];
    end
end

% Create UI figure for piecewise function input
piece_fig = uifigure('Name', [disturbance_label ' Piecewise Function Input'], 'Position', [400, 200, 800, 600]);

% Add plot area at top
ax = uiaxes(piece_fig, 'Position', [150, 350, 500, 230]);
title(ax, 'Piecewise Function Preview');
xlabel(ax, 'Time (s)');
ylabel(ax, disturbance_label);
grid(ax, 'on');

% Create table for sections with input fields
columnNames = {'Start Time', 'End Time', 'Time Steps', disturbance_label};
columnTypes = {'numeric', 'numeric', 'numeric', 'char'};
columnEditable = [true, true, true, true];
columnWidth = {80, 80, 80, 350};

% Initialize data with default values
data = cell(num_sections, 4);
for i = 1:num_sections
    if i == 1
        data{i, 1} = 0; % Default start time for first section
    else
        data{i, 1} = data{i-1, 2}; % Start time = end time of previous section
    end
    data{i, 2} = data{i, 1} + 10; % Default end time
    data{i, 3} = 0.1; % Default time steps
    data{i, 4} = ''; % Empty function by default
end

% Create the table
pieceTable = uitable(piece_fig, 'Position', [150, 200, 500, 120], ...
    'Data', data, 'ColumnName', columnNames, 'ColumnEditable', columnEditable, ...
    'ColumnWidth', columnWidth);

% Function info text
uilabel(piece_fig, 'Text', 'Enter any valid MATLAB expression using "t" as the variable.', ...
    'Position', [150, 170, 500, 22]);
uilabel(piece_fig, 'Text', 'Examples: sin(t), t^2, sqrt(t), pi*cos(t), exp(-t)', ...
    'Position', [150, 150, 500, 22]);
uilabel(piece_fig, 'Text', 'Note: For section boundaries, make sure end time of section i equals start time of section i+1', ...
    'Position', [150, 130, 500, 22]);

% Buttons at bottom
buttonWidth = 100;
buttonHeight = 30;
panelWidth = piece_fig.Position(3);
buttonY = 60;

% Preview button
previewButton = uibutton(piece_fig, 'Position', [(panelWidth/2)-150, buttonY, buttonWidth, buttonHeight], ...
    'Text', 'Preview', ...
    'ButtonPushedFcn', @(btn, event) disturbance_preview_piecewise_function(ax, pieceTable, disturbance_label));

% Confirm button
confirmButton = uibutton(piece_fig, 'Position', [(panelWidth/2)-50, buttonY, buttonWidth, buttonHeight], ...
    'Text', 'Confirm', ...
    'ButtonPushedFcn', @(btn, event) disturbance_confirm_piecewise_function(selection_fig, piece_fig, pieceTable));

% Cancel button
cancelButton = uibutton(piece_fig, 'Position', [(panelWidth/2)+50, buttonY, buttonWidth, buttonHeight], ...
    'Text', 'Cancel', ...
    'ButtonPushedFcn', @(btn, event) disturbance_cancel_input(piece_fig));

% Handle window close event
piece_fig.CloseRequestFcn = @(src, event) disturbance_cancel_input(src);

% Wait for user to confirm or cancel
uiwait(piece_fig);
end