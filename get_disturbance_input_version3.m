function [disturbance, num_sections, start_time, end_time] = get_disturbance_input_version3(disturbance_label)
% GET_DISTURBANCE_INPUT_VERSION3 - Get disturbance value with enhanced UI
% Similar to get_target_value_version2 but for disturbance inputs
%
% Parameters:
%   disturbance_label - String label for the disturbance (e.g., 'd_1' or 'd_2')
%
% Returns:
%   disturbance: Function handle (continuous) or cell array of function handles (piecewise)
%   num_sections: Number of sections (1 for continuous)
%   start_time: Vector of section start times
%   end_time: Vector of section end times

% Initialize return values
disturbance = [];
num_sections = [];
start_time = [];
end_time = [];

% Create a figure for user interaction
selection_fig = uifigure('Name', ['Select ' disturbance_label ' Function Type'], 'Position', [500, 500, 300, 150]);

% Set default values in case window is closed
setappdata(selection_fig, 'disturbance', []);
setappdata(selection_fig, 'num_sections', []);
setappdata(selection_fig, 'start_time', []);
setappdata(selection_fig, 'end_time', []);
setappdata(selection_fig, 'time_steps', []);

% Add buttons for continuous and piecewise functions
uibutton(selection_fig, 'Position', [50, 80, 200, 50], 'Text', 'Continuous Function', ...
    'ButtonPushedFcn', @(btn, event) disturbance_continuous_function_input(selection_fig, disturbance_label));
uibutton(selection_fig, 'Position', [50, 20, 200, 50], 'Text', 'Define Function in Sections', ...
    'ButtonPushedFcn', @(btn, event) disturbance_piecewise_function_input(selection_fig, disturbance_label));

% Handle window close event
selection_fig.CloseRequestFcn = @(src, event) disturbance_selection_figure_closed(src);

% Wait for user input
uiwait(selection_fig);

% Check if figure still exists
if isvalid(selection_fig)
    % Get values from app data
    disturbance = getappdata(selection_fig, 'disturbance');
    num_sections = getappdata(selection_fig, 'num_sections');
    start_time = getappdata(selection_fig, 'start_time');
    end_time = getappdata(selection_fig, 'end_time');
    time_steps = getappdata(selection_fig, 'time_steps');
    
    % If values were successfully obtained from UI, pass them to get_time_vector
    if ~isempty(start_time) && ~isempty(end_time) && ~isempty(time_steps) && ~isempty(num_sections)
        % Store the time vector parameters with a unique name for this disturbance
        % This uses a modified version of get_time_vector that doesn't share persistent variables
        [start_time, end_time, time_steps] = get_time_vector_disturbance(num_sections, start_time, end_time, time_steps, disturbance_label);
    end
    
    % Close the selection figure
    delete(selection_fig);
end
end

% The helper functions are moved to separate files