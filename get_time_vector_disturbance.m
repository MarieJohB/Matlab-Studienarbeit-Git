function [start_time, end_time, time_steps] = get_time_vector_disturbance(num_sections, varargin)
% GET_TIME_VECTOR_DISTURBANCE - Get time vector parameters for disturbance inputs
% Similar to get_time_vector but uses separate persistent variables for each disturbance
%
% Parameters:
%   num_sections - Number of time sections
%   varargin - Optional parameters: 
%     If 3 additional args: start_time, end_time, time_steps
%     If 4 additional args: start_time, end_time, time_steps, disturbance_label
% Returns:
%   start_time - Vector of start times for each section
%   end_time - Vector of end times for each section
%   time_steps - Vector of time steps for each section

% Determine the disturbance label (default to 'd_generic' if not provided)
disturbance_label = 'd_generic';
if nargin >= 5
    disturbance_label = varargin{4};
end

% Use persistent containers.Map objects to store values for each disturbance
persistent start_values end_values steps_values;
if isempty(start_values)
    start_values = containers.Map();
    end_values = containers.Map();
    steps_values = containers.Map();
end

% Check if parameters were passed directly
if nargin >= 4
    % Use parameters provided directly
    start_time = varargin{1};
    end_time = varargin{2};
    time_steps = varargin{3};
    
    % Store these values for this disturbance
    start_values(disturbance_label) = start_time;
    end_values(disturbance_label) = end_time;
    steps_values(disturbance_label) = time_steps;
    
    return;
end

% Check if we already have values for this disturbance
if isKey(start_values, disturbance_label) && isKey(end_values, disturbance_label) && ...
   isKey(steps_values, disturbance_label) && length(start_values(disturbance_label)) == num_sections
    
    % Use the stored values
    start_time = start_values(disturbance_label);
    end_time = end_values(disturbance_label);
    time_steps = steps_values(disturbance_label);
else
    % If no stored values or section count changed, prompt for new ones
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
    
    % Store values for this disturbance
    start_values(disturbance_label) = start_time;
    end_values(disturbance_label) = end_time;
    steps_values(disturbance_label) = time_steps;
end
end

function [valid] = validate_time_parameters_disturbance(start_time, end_time, time_steps)
% VALIDATE_TIME_PARAMETERS_DISTURBANCE - Validate time parameters for section creation
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