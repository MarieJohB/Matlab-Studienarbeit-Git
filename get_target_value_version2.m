
function target_value = get_target_value_version2()
    % Initialize target_value to be empty
    target_value = [];

    % Create a figure for user interaction
    selection_fig = uifigure('Name', 'Select Function Type', 'Position', [500, 500, 300, 150]);
    
    % Add buttons for continuous and piecewise functions
    uibutton(selection_fig, 'Position', [50, 80, 200, 50], 'Text', 'Continuous Function', ...
        'ButtonPushedFcn', @(btn, event) continuous_function_input(selection_fig));
    uibutton(selection_fig, 'Position', [50, 20, 200, 50], 'Text', 'Define Function in Sections', ...
        'ButtonPushedFcn', @(btn, event) piecewise_function_input(selection_fig));
end


function continuous_function_input(selection_fig)

    % Close the selection figure
    close(selection_fig);

    % Initialize target_value to be empty
    target_value = [];
      % Loop until valid input is provided or user cancels the dialog
    while isempty(target_value)
        % Prompt for target value (constant or time-dependent)
        prompt = {'Enter the target value (use "t" for time-dependent values, e.g., "5" or "sin(t)"): '};
        dlgtitle = 'Target Value Input';
        dims = [1 50];
        answer = inputdlg(prompt, dlgtitle, dims);
        % Check if user cancelled the dialog
        if isempty(answer)
            disp('Operation cancelled by user.');
            return;
        end
        % Validate input
        try
            % Replace commas with periods
            target_value_str = strrep(answer{1}, ',', '.');
            % Test the input
            t = 0; % Dummy value for validation
            target_value = eval(target_value_str);
            if isnan(target_value)
                error('Invalid input. Please enter a valid expression.');
            end
            % Define the target value as a function of time if it contains "t"
            if contains(target_value_str, 't')
                target_value = str2func(['@(t)', target_value_str]);
            else
                % If it's a constant value, create a function that returns the constant + 0*t
                target_value = str2func(['@(t)', target_value_str, ' + 0*t']);
            end
        catch
            uiwait(msgbox('Invalid input. Please enter the expression correctly.', 'Error', 'error'));
            target_value = [];
        end
    end
    % Display the target value type
    if ~isempty(target_value)
        disp('The target value has been set.');
    end
    
end

function piecewise_function_input(selection_fig)
    % Close the selection figure
    close(selection_fig);
    
    % initalize variable for number of sections
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
    
    % Get the number of sections
    % maximum number of sections is 5
    % number cannot be decimal

    % Replace commas with periods
    num_sections_str = strrep(answer{1}, ',', '.');
    num_sections = str2double(num_sections_str);

    if isnan(num_sections) || num_sections < 1 || num_sections > 5 || mod(num_sections, 1) ~= 0
        uiwait(msgbox('Invalid input. Please enter a number between 1 and 5.', 'Error', 'error'));
        num_sections = [];
    end

    end



        % initalize variables to for  start time, end time and steps
        start_time = NaN * ones(1, num_sections);
        end_time = NaN * ones(1, num_sections);
        time_steps = NaN * ones(1, num_sections);
        target_value = NaN * ones(1, num_sections);

    for i = 1:num_sections
    
        while isnan(start_time(i)) || isnan(end_time(i)) || isnan(time_steps(i)) || isnan(target_value(i))
            
            if i==1
             
            % define request for user input
            prompt = {'Start Time:', 'End Time:', 'Time Steps:', 'Function for this section:'};
            dlgtitle = 'Input';
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

            %checking if end > start:
            if start_time(i) > end_time(i)
                uiwait(msgbox('Invalid input. Please make sure that end happens after start.', 'Error', 'error'));
                start_time(i) = NaN;
                end_time(i) = NaN;
            end

            % checking if step size is compatible with time span
            numSteps = (end_time(i) - start_time(i)) / time_steps(i);
            if mod(numSteps, 1) ~= 0
                uiwait(msgbox('Invalid input. Please enter suitable size for time steps.', 'Error', 'error'));
                time_steps(i) = NaN;
            end


            % Validate input
            try
            % Replace commas with periods
            target_value_str = strrep(answer{1}, ',', '.');
            % Test the input
            t = 0; % Dummy value for validation
            target_value(i) = eval(target_value_str);
            if isnan(target_value(i))
                error('Invalid input. Please enter a valid expression.');
            end
            % Define the target value as a function of time if it contains "t"
            if contains(target_value_str, 't')
                target_value(i) = str2func(['@(t)', target_value_str]);
            else
                % If it's a constant value, create a function that returns the constant + 0*t
                target_value(i) = str2func(['@(t)', target_value_str, ' + 0*t']);
            end
            catch
                uiwait(msgbox('Invalid input. Please enter the expression correctly.', 'Error', 'error'));
                target_value(i) = NaN;
            end
        




            elseif i > 1 
                % for the following sections
                % end time is the start time for the next section
                % define request for user input
                prompt = {'End Time:', 'Time Steps:', 'Function for this section:'};
                dlgtitle = 'Input';
                dims = [1 35];
                answer = inputdlg(prompt, dlgtitle, dims);

                % Check if user cancelled the dialog
                if isempty(answer)
                    disp('Operation cancelled by user.');
                       return;
                end

            
                start_time(i) = end_time(i-1);

                end_time_str = strrep(answer{1}, ',', '.');
                end_time(i) = str2double(end_time_str);
    
                time_steps_str = strrep(answer{2}, ',', '.');
                time_steps(i) = str2double(time_steps_str);

                %checking if end > start:
                if start_time(i) > end_time(i)
                    uiwait(msgbox('Invalid input. Please make sure that end happens after start.', 'Error', 'error'));
                    start_time(i) = [];
                    end_time(i) = [];
                end

                % checking if step size is compatible with time span
                numSteps = (end_time(i) - start_time(i)) / time_steps(i);
                if mod(numSteps, 1) ~= 0
                    uiwait(msgbox('Invalid input. Please enter suitable size for time steps.', 'Error', 'error'));
                    time_steps(i) = [];
                end


                % Validate input
                try
                % Replace commas with periods
                target_value_str = strrep(answer{3}, ',', '.');
                % Test the input
                t = 0; % Dummy value for validation
                target_value(i) = eval(target_value_str);
                if isnan(target_value(i))
                error('Invalid input. Please enter a valid expression.');
                end
                % Define the target value as a function of time if it contains "t"
                if contains(target_value_str, 't')
                    target_value(i) = str2func(['@(t)', target_value_str]);
                    else
                    % If it's a constant value, create a function that returns the constant + 0*t
                    target_value(i) = str2func(['@(t)', target_value_str, ' + 0*t']);
                end
                catch
                uiwait(msgbox('Invalid input. Please enter the expression correctly.', 'Error', 'error'));
                target_value(i) = NaN;
                end
        
            end
        
        end
    end

disp('Start times:'); 
disp(start_time);
disp('End times:'); 
disp(end_time);
disp('Time steps:'); 
disp(time_steps);
disp('Target values:'); 
disp(target_value);

end





