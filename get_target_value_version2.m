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

    % Create a figure for continuous function input
    function_fig = uifigure('Name', 'Continuous Target Value Input', 'Position', [500, 500, 300, 150]);
    
    prompt = {'Enter the continuous target value (e.g., "5" or "sin(t)"): '};
    dlgtitle = 'Continuous Target Value Input';
    dims = [1 50];
    answer = inputdlg(prompt, dlgtitle, dims);
    
    % Check if user cancelled the dialog
    if isempty(answer)
        disp('Operation cancelled by user.');
        return;
    end
    
    % Replace commas with dots
    target_value_str = strrep(answer{1}, ',', '.');
    
    % Test the input
    try
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
        
        disp('The continuous target value has been set.');
    catch
        uiwait(msgbox('Invalid input. Please enter the expression correctly.', 'Error', 'error'));
        target_value = [];
    end

    % Close the function figure
    close(function_fig);
end

function piecewise_function_input(selection_fig)
    % Close the selection figure
    close(selection_fig);
    
    % Create a figure for selecting number of sections
    fig = uifigure('Name', 'Select Number of Sections', 'Position', [500, 500, 300, 150]);
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
    num_sections = str2double(answer{1});
    if isnan(num_sections) || num_sections < 1 || num_sections > 5
        uiwait(msgbox('Invalid input. Please enter a number between 1 and 5.', 'Error', 'error'));
        return;
    end
    
    % Initialize section inputs
    section_inputs = cell(num_sections, 1);
    
    % Loop through sections and get input for each section
    for i = 1:num_sections
        fig = uifigure('Name', ['Section ', num2str(i)], 'Position', [500, 500, 300, 250]);
        uilabel(fig, 'Position', [20, 180, 100, 22], 'Text', ['Section ', num2str(i), ':']);
        
        if i == 1
            start_time_label = uilabel(fig, 'Position', [20, 140, 100, 22], 'Text', 'Start time:');
            start_time_edit = uieditfield(fig, 'numeric', 'Position', [120, 140, 150, 22]);
        else
            start_time_label = uilabel(fig, 'Position', [20, 140, 100, 22], 'Text', ['Start time (t', num2str(i-1), '):']);
            start_time_edit = uieditfield(fig, 'numeric', 'Position', [120, 140, 150, 22], 'Value', section_inputs{i-1}.end_time);
        end
        
        end_time_label = uilabel(fig, 'Position', [20, 100, 100, 22], 'Text', 'End time:');
        end_time_edit = uieditfield(fig, 'numeric', 'Position', [120, 100, 150, 22]);
        
        value_label = uilabel(fig, 'Position', [20, 60, 100, 22], 'Text', 'Value:');
        value_edit = uieditfield(fig, 'text', 'Position', [120, 60, 150, 22]);
        
        uibutton(fig, 'Position', [100, 20, 100, 30], 'Text', 'Confirm', ...
            'ButtonPushedFcn', @(btn, event) confirm_section(fig, i, start_time_edit.Value, end_time_edit.Value, value_edit.Value, section_inputs));
        
        uiwait(fig);
    end
    
    % Combine all sections into one piecewise function
    sections = "";
    for i = 1:num_sections
        if i > 1
            sections = sections + " ";
        end
        sections = sections + "t >= " + num2str(section_inputs{i}.start_time) + ", " + section_inputs{i}.value + ", t < " + num2str(section_inputs{i}.end_time);
    end
    target_value_str = "piecewise(" + sections + ")";
    target_value = eval(['@(t) ' target_value_str]);
    disp('The piecewise target value has been set.');
end

function confirm_section(fig, section_number, start_time, end_time, value, section_inputs)
    % Validate the inputs
    if isempty(start_time) || isempty(end_time) || isempty(value) || start_time >= end_time
        uiwait(msgbox('Invalid input. Please ensure all fields are filled and start time is less than end time.', 'Error', 'error'));
        return;
    end
    
    % Save the inputs
    section_inputs{section_number} = struct('start_time', start_time, 'end_time', end_time, 'value', value);
    
    % Close the figure
    close(fig);
end

