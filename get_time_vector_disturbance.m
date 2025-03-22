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
        % Use parameters provided from UI
        start_time = varargin{1};
        end_time = varargin{2};
        time_steps = varargin{3};
        
        % Store these values in the map for this disturbance
        start_values(disturbance_label) = start_time;
        end_values(disturbance_label) = end_time;
        steps_values(disturbance_label) = time_steps;
        
        return;
    end

    % Check if we already have values for this disturbance
    if isKey(start_values, disturbance_label) && isKey(end_values, disturbance_label) && isKey(steps_values, disturbance_label)
        % Use the stored values
        start_time = start_values(disturbance_label);
        end_time = end_values(disturbance_label);
        time_steps = steps_values(disturbance_label);
    else
        % If no stored values, prompt the user for new ones
        % Initialize variables for start time, end time and steps
        start_time = NaN * ones(1, num_sections);
        end_time = NaN * ones(1, num_sections);
        time_steps = NaN * ones(1, num_sections);
        
        for i = 1:num_sections
            while isnan(start_time(i)) || isnan(end_time(i)) || isnan(time_steps(i))
                if i==1
                    % Define request for user input
                    prompt = {'Start Time:', 'End Time:', 'Time Steps:'};
                    dlgtitle = [disturbance_label ' Time Input'];
                    dims = [1 35];
                    answer = inputdlg(prompt, dlgtitle, dims);
                    
                    % Check if user cancelled the dialog
                    if isempty(answer)
                        disp('Operation cancelled by user.');
                        return;
                    end
                    
                    start_time_str = strrep(answer{1}, ',', '.');
                    start_time(i) = str2double(start_time_str);
                    
                    end_time_str = strrep(answer{2}, ',', '.');
                    end_time(i) = str2double(end_time_str);
                    
                    time_steps_str = strrep(answer{3}, ',', '.');
                    time_steps(i) = str2double(time_steps_str);
                    
                    % Checking if end > start
                    if start_time(i) >= end_time(i)
                        uiwait(msgbox('Invalid input. Please make sure that end happens after start.', 'Error', 'error'));
                        start_time(i) = NaN;
                        end_time(i) = NaN;
                    end
                    
                    % Check if end and start > 0
                    if start_time(i) < 0
                        uiwait(msgbox('Invalid input. Please enter positive values only.', 'Error', 'error'));
                        start_time(i) = NaN;
                    end
                    
                    % Checking if step size is compatible with time span
                    numSteps = (end_time(i) - start_time(i)) / time_steps(i);
                    if mod(numSteps, 1) ~= 0
                        uiwait(msgbox('Invalid input. Please enter suitable size for time steps.', 'Error', 'error'));
                        time_steps(i) = NaN;
                    end
                    
                elseif i > 1 
                    % For the following sections
                    % End time is the start time for the next section
                    start_time(i) = end_time(i-1);
                    
                    % Adding information for user: start time of current section
                    message = sprintf('Start time of this section = %f. \nPlease enter end time and steps for this section', start_time(i));
                    uiwait(msgbox(message, 'Information', 'modal'));
                    
                    % Define request for user input
                    prompt = {'End Time:', 'Time Steps:'};
                    dlgtitle = [disturbance_label ' Section ' num2str(i) ' Input'];
                    dims = [1 35];
                    
                    answer = inputdlg(prompt, dlgtitle, dims);
                    
                    % Check if user cancelled the dialog
                    if isempty(answer)
                        disp('Operation cancelled by user.');
                        return;
                    end
                    
                    end_time_str = strrep(answer{1}, ',', '.');
                    end_time(i) = str2double(end_time_str);
                    
                    time_steps_str = strrep(answer{2}, ',', '.');
                    time_steps(i) = str2double(time_steps_str);
                    
                    % Checking if end > start
                    if start_time(i) > end_time(i)
                        uiwait(msgbox('Invalid input. Please make sure that end happens after start.', 'Error', 'error'));
                        start_time(i) = NaN;
                        end_time(i) = NaN;
                    end
                    
                    % Checking if step size is compatible with time span
                    numSteps = (end_time(i) - start_time(i)) / time_steps(i);
                    if mod(numSteps, 1) ~= 0
                        uiwait(msgbox('Invalid input. Please enter suitable size for time steps.', 'Error', 'error'));
                        time_steps(i) = NaN;
                    end
                end
            end
        end
        
        % Save the defined values in persistent variables for this disturbance
        start_values(disturbance_label) = start_time;
        end_values(disturbance_label) = end_time;
        steps_values(disturbance_label) = time_steps;
    end
end