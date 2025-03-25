function [target_value, num_sections, start_time, end_time] = get_target_value_version2()
% GET_TARGET_VALUE_VERSION2 - Enhanced UI for getting target value for the reference signal

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

% Initialize return values
target_value = [];
num_sections = [];
start_time = [];
end_time = [];

% Create a selection figure with improved styling
selection_fig = uifigure('Name', 'Reference Signal Type', 'Position', [500, 300, 400, 300]);
selection_fig.Color = appColors.background;

% Add title panel with enhanced styling
titlePanel = uipanel(selection_fig, 'Position', [10 250 380 40], 'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
titleLabel = uilabel(titlePanel, 'Text', 'Select Reference Signal Type', ...
    'Position', [0 0 380 40], 'FontSize', 16, 'FontWeight', 'bold', ...
    'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');

% Set default values in case window is closed
setappdata(selection_fig, 'target_value', []);
setappdata(selection_fig, 'num_sections', []);
setappdata(selection_fig, 'start_time', []);
setappdata(selection_fig, 'end_time', []);
setappdata(selection_fig, 'time_steps', []);

% Buttons directly in figure (no panel)
contBtn = uibutton(selection_fig, 'push', 'Position', [100, 180, 200, 50], ...
    'Text', 'Continuous Function', ...
    'FontSize', 12, 'FontWeight', 'bold', ...
    'BackgroundColor', appColors.buttonPrimary, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) target_continuous_function_input(selection_fig));

pieceBtn = uibutton(selection_fig, 'push', 'Position', [100, 120, 200, 50], ...
    'Text', 'Piecewise Function', ...
    'FontSize', 12, 'FontWeight', 'bold', ...
    'BackgroundColor', appColors.buttonPrimary, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) target_piecewise_function_input(selection_fig));

cancelBtn = uibutton(selection_fig, 'push', 'Position', [150, 60, 100, 30], ...
    'Text', 'Cancel', ...
    'FontSize', 12, ...
    'BackgroundColor', appColors.buttonCancel, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) closeSelectionFig());

% Handle window close event
selection_fig.CloseRequestFcn = @(src, event) closeSelectionFig();

% Function to close the selection figure
function closeSelectionFig()
    uiresume(selection_fig);
    delete(selection_fig);
end

% Wait for user input
uiwait(selection_fig);

% Check if figure still exists
if isvalid(selection_fig)
    % Get values from app data
    target_value = getappdata(selection_fig, 'target_value');
    num_sections = getappdata(selection_fig, 'num_sections');
    start_time = getappdata(selection_fig, 'start_time');
    end_time = getappdata(selection_fig, 'end_time');
    time_steps = getappdata(selection_fig, 'time_steps');
    
    % If values were successfully obtained from UI, pass them to get_time_vector
    if ~isempty(start_time) && ~isempty(end_time) && ~isempty(time_steps) && ~isempty(num_sections)
        % Call get_time_vector with parameters to store the values
        [start_time, end_time, time_steps] = get_time_vector(num_sections, start_time, end_time, time_steps);
    end
    
    % Close the selection figure
    delete(selection_fig);
end
end