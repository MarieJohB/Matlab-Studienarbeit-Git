function [target_value, num_sections, start_time, end_time] = get_target_value_version2()
% GET_TARGET_VALUE_VERSION2 - Get target value for the reference signal
% This function creates a main selection dialog where users can choose between
% continuous or piecewise functions, then calls the appropriate helper function.
%
% Returns:
%   target_value: Function handle (continuous) or cell array of function handles (piecewise)
%   num_sections: Number of sections (1 for continuous)
%   start_time: Vector of section start times
%   end_time: Vector of section end times

% Initialize return values
target_value = [];
num_sections = [];
start_time = [];
end_time = [];

% Create a selection figure with improved styling
selection_fig = uifigure('Name', 'Select Function Type', 'Position', [500, 500, 300, 150]);

% Set default values in case window is closed
setappdata(selection_fig, 'target_value', []);
setappdata(selection_fig, 'num_sections', []);
setappdata(selection_fig, 'start_time', []);
setappdata(selection_fig, 'end_time', []);
setappdata(selection_fig, 'time_steps', []);

% Add buttons for continuous and piecewise functions
uibutton(selection_fig, 'Position', [50, 80, 200, 50], 'Text', 'Continuous Function', ...
    'ButtonPushedFcn', @(btn, event) target_continuous_function_input(selection_fig));
uibutton(selection_fig, 'Position', [50, 20, 200, 50], 'Text', 'Define Function in Sections', ...
    'ButtonPushedFcn', @(btn, event) target_piecewise_function_input(selection_fig));

% Handle window close event
selection_fig.CloseRequestFcn = @(src, event) target_selection_figure_closed(src);

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

% The helper functions are moved to separate files