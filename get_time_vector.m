function [start_time, end_time, time_steps] = get_time_vector(num_sections, varargin)
% GET_TIME_VECTOR - Get time vector parameters with optional direct input and better integration
% 
% Parameters:
%   num_sections - Number of time sections
%   varargin - Optional parameters: start_time, end_time, time_steps
% 
% Returns:
%   start_time - Vector of start times for each section
%   end_time - Vector of end times for each section
%   time_steps - Vector of time steps for each section

% Initialize persistent variables to cache time vector parameters
persistent start_value;
persistent end_value;
persistent steps_value;

% Check if parameters were passed directly (indicates values from UI)
if nargin >= 4
    % Use parameters provided from UI
    start_value = varargin{1};
    end_value = varargin{2};
    steps_value = varargin{3};
    
    % Return the values directly
    start_time = start_value;
    end_time = end_value;
    time_steps = steps_value;
    return;
end

% If no stored values, or parameters have been cleared, create new ones
if isempty(start_value) || isempty(end_value) || isempty(steps_value) || ...
   length(start_value) ~= num_sections || length(end_value) ~= num_sections
    
    % Initialize vectors for time parameters
    start_time = zeros(1, num_sections);
    end_time = zeros(1, num_sections);
    time_steps = zeros(1, num_sections);
    
    % Get input for each section
    for i = 1:num_sections
        if i == 1
            % First section: get all parameters
            prompt = {'Start Time:', 'End Time:', 'Time Steps:'};
            dlgtitle = 'Time Parameters';
            dims = [1 35];
            defaultans = {'0', '10', '0.1'};
            answer = inputdlg(prompt, dlgtitle, dims, defaultans);
            
            % Check if user cancelled
            if isempty(answer)
                disp('Operation cancelled by user.');
                start_time = [];
                end_time = [];
                time_steps = [];
                return;
            end
            
            % Process input (allowing decimal commas)
            start_time_str = strrep(answer{1}, ',', '.');
            end_time_str = strrep(answer{2}, ',', '.');
            time_steps_str = strrep(answer{3}, ',', '.');
            
            start_time(i) = str2double(start_time_str);
            end_time(i) = str2double(end_time_str);
            time_steps(i) = str2double(time_steps_str);
            
        else
            % Subsequent sections: start time = end time of previous section
            start_time(i) = end_time(i-1);
            
            % Show current section info
            message = sprintf('Start time of section %d = %.4f', i, start_time(i));
            
            % Get end time and time steps
            prompt = {'End Time:', 'Time Steps:'};
            dlgtitle = sprintf('Section %d Parameters', i);
            dims = [1 35];
            defaultans = {num2str(start_time(i) + 10), '0.1'};
            answer = inputdlg(prompt, dlgtitle, dims, defaultans, 'on');
            
            % Check if user cancelled
            if isempty(answer)
                disp('Operation cancelled by user.');
                start_time = [];
                end_time = [];
                time_steps = [];
                return;
            end
            
            % Process input (allowing decimal commas)
            end_time_str = strrep(answer{1}, ',', '.');
            time_steps_str = strrep(answer{2}, ',', '.');
            
            end_time(i) = str2double(end_time_str);
            time_steps(i) = str2double(time_steps_str);
        end
        
        % Validate input for this section
        if ~validateTimeParameters(start_time(i), end_time(i), time_steps(i))
            % Reset and try again if validation fails
            i = i - 1;
            continue;
        end
    end
    
    % Store values in persistent variables
    start_value = start_time;
    end_value = end_time;
    steps_value = time_steps;
else
    % Return previously stored values
    start_time = start_value;
    end_time = end_value;
    time_steps = steps_value;
end
end

function [t] = create_linear_time_vector(num_sections)
% CREATE_LINEAR_TIME_VECTOR - Creates a linear time vector for all sections
%
% Parameters:
%   num_sections - Number of time sections
% Returns:
%   t - Combined linear time vector for all sections

% Initialize time vector
t = [];

% Check if num_sections is valid
if isempty(num_sections) || num_sections < 1
    disp('Invalid number of sections.');
    return;
end

% Get time parameters
[start_time, end_time, time_steps] = get_time_vector(num_sections);

% Check if time parameters were returned
if isempty(start_time) || isempty(end_time) || isempty(time_steps)
    disp('Operation cancelled: Time parameters not provided.');
    return;
end

% Create time vector for each section
for i = 1:num_sections
    % Calculate number of points in this section
    num_points = round((end_time(i) - start_time(i)) / time_steps(i)) + 1;
    
    % Create time vector for this section
    section_t = linspace(start_time(i), end_time(i), num_points);
    
    % If this is not the first section, remove the first point to avoid duplication
    if i > 1 && ~isempty(t) && ~isempty(section_t)
        section_t = section_t(2:end);
    end
    
    % Append to overall time vector
    t = [t, section_t];
end

% Display time vector information
if ~isempty(t)
    disp(['Created time vector with ', num2str(length(t)), ' points from ', ...
          num2str(t(1)), ' to ', num2str(t(end)), ' across ', num2str(num_sections), ' section(s).']);
else
    disp('Warning: Created an empty time vector.');
end
end

function [values_vector] = evaluate_piecewise_function(function_vector, time_vector, start_times, end_times)
% EVALUATE_PIECEWISE_FUNCTION - Evaluates a piecewise function for a given time vector
%
% Parameters:
%   function_vector - Function handle (continuous) or cell array of function handles (piecewise)
%   time_vector - Time vector for evaluation
%   start_times - Vector of section start times
%   end_times - Vector of section end times
%
% Returns:
%   values_vector - Evaluated function values

% Initialize values vector
values_vector = zeros(size(time_vector));

% Check if continuous or piecewise function
if ~iscell(function_vector)
    % Continuous function case
    if isa(function_vector, 'function_handle')
        values_vector = arrayfun(function_vector, time_vector);
    else
        % Constant value
        values_vector = function_vector * ones(size(time_vector));
    end
else
    % Piecewise function case
    num_sections = length(function_vector);
    
    for i = 1:length(time_vector)
        current_time = time_vector(i);
        
        % Find the active section for this time point
        section_idx = find(current_time >= start_times & current_time <= end_times, 1, 'first');
        
        % Handle points outside defined sections (use nearest section)
        if isempty(section_idx)
            if current_time < start_times(1)
                section_idx = 1;
            elseif current_time > end_times(end)
                section_idx = num_sections;
            end
        end
        
        % Evaluate function for current time point
        if ~isempty(section_idx) && section_idx <= num_sections
            if isa(function_vector{section_idx}, 'function_handle')
                values_vector(i) = function_vector{section_idx}(current_time);
            else
                values_vector(i) = function_vector{section_idx};
            end
        end
    end
end
end

function [valid] = validateTimeParameters(start_time, end_time, time_steps)
% VALIDATETIMEPARAMETERS - Validates time parameters for section creation
%
% Parameters:
%   start_time - Section start time
%   end_time - Section end time
%   time_steps - Section time step size
%
% Returns:
%   valid - True if parameters are valid, false otherwise

valid = true;

% Check if end time is greater than start time
if end_time <= start_time
    uiwait(msgbox('End time must be greater than start time.', 'Error', 'error'));
    valid = false;
    return;
end

% Check if start time is non-negative
if start_time < 0
    uiwait(msgbox('Start time must be non-negative.', 'Error', 'error'));
    valid = false;
    return;
end

% Check if time step is positive
if time_steps <= 0
    uiwait(msgbox('Time step must be positive.', 'Error', 'error'));
    valid = false;
    return;
end

% Check if time step divides the interval evenly
numSteps = (end_time - start_time) / time_steps;
if abs(round(numSteps) - numSteps) > 1e-6
    uiwait(msgbox('Time step must divide the interval evenly.', 'Error', 'error'));
    valid = false;
    return;
end

end