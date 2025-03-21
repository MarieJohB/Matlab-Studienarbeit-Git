function [target_value, num_sections, start_time, end_time] = get_target_value_version2()
    % GET_TARGET_VALUE_VERSION2 - Get target value with enhanced UI
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

    % Create a figure for user interaction
    selection_fig = uifigure('Name', 'Select Function Type', 'Position', [500, 500, 300, 150]);
    
    % Set default values in case window is closed
    setappdata(selection_fig, 'target_value', []);
    setappdata(selection_fig, 'num_sections', []);
    setappdata(selection_fig, 'start_time', []);
    setappdata(selection_fig, 'end_time', []);
    setappdata(selection_fig, 'time_steps', []);
    
    % Add buttons for continuous and piecewise functions
    uibutton(selection_fig, 'Position', [50, 80, 200, 50], 'Text', 'Continuous Function', ...
        'ButtonPushedFcn', @(btn, event) continuous_function_input(selection_fig));
    uibutton(selection_fig, 'Position', [50, 20, 200, 50], 'Text', 'Define Function in Sections', ...
        'ButtonPushedFcn', @(btn, event) piecewise_function_input(selection_fig));

    % Handle window close event
    selection_fig.CloseRequestFcn = @(src, event) selectionFigureClosed(src);
    
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

function selectionFigureClosed(fig)
    % Handle the case when the figure is closed without selection
    uiresume(fig);
end

function continuous_function_input(selection_fig)
    % Continuous function input with integrated time vector input
    
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
    
    % Time Steps - now to the right of End Time
    uilabel(pnl, 'Text', 'Time Steps:', 'Position', [255, 5, 60, 22]);
    stepField = uieditfield(pnl, 'numeric', 'Position', [320, 5, 50, 22], 'Value', 0.1);
    
    % Function input - r(t) label left-aligned to input field
    uilabel(cont_fig, 'Text', 'r(t):', 'Position', [150, 100, 40, 22]);
    funcField = uieditfield(cont_fig, 'text', 'Position', [190, 100, 260, 22], 'Value', '');
    
    % Function info text
    uilabel(cont_fig, 'Text', 'Enter any valid MATLAB expression using "t" as the variable.', ...
        'Position', [100, 70, 400, 22]);
    uilabel(cont_fig, 'Text', 'Examples: sin(t), t^2, sqrt(t), pi*cos(t), exp(-t)', ...
        'Position', [100, 50, 400, 22]);
    
    % Buttons centered at bottom
    buttonWidth = 100;
    buttonHeight = 25;
    panelWidth = cont_fig.Position(3);
    buttonY = 15;
    
    % Preview button
    previewButton = uibutton(cont_fig, 'Position', [(panelWidth/2)-150, buttonY, buttonWidth, buttonHeight], ...
        'Text', 'Preview', ...
        'ButtonPushedFcn', @(btn, event) previewContinuousFunction(ax, startField, endField, stepField, funcField));
    
    % Confirm button
    confirmButton = uibutton(cont_fig, 'Position', [(panelWidth/2)-50, buttonY, buttonWidth, buttonHeight], ...
        'Text', 'Confirm', ...
        'ButtonPushedFcn', @(btn, event) confirmContinuousFunction(selection_fig, cont_fig, startField, endField, stepField, funcField));
    
    % Cancel button
    cancelButton = uibutton(cont_fig, 'Position', [(panelWidth/2)+50, buttonY, buttonWidth, buttonHeight], ...
        'Text', 'Cancel', ...
        'ButtonPushedFcn', @(btn, event) cancelInput(cont_fig));
    
    % Handle window close event
    cont_fig.CloseRequestFcn = @(src, event) cancelInput(src);
    
    % Wait for user to confirm or cancel
    uiwait(cont_fig);
end

function piecewise_function_input(selection_fig)
    % Piecewise function input with integrated time vector input
    
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
    piece_fig = uifigure('Name', 'Piecewise Function Input', 'Position', [400, 200, 800, 600]);
    
    % Add plot area at top, centered
    ax = uiaxes(piece_fig, 'Position', [150, 350, 500, 230]);
    title(ax, 'Piecewise Function Preview');
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'r(t)');
    grid(ax, 'on');
    
    % Create table for sections with input fields
    columnNames = {'Start Time', 'End Time', 'Time Steps', 'r(t)'};
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
    
    % Create the table - centered below plot
    pieceTable = uitable(piece_fig, 'Position', [150, 200, 500, 120], ...
        'Data', data, 'ColumnName', columnNames, 'ColumnEditable', columnEditable, ...
        'ColumnWidth', columnWidth);
    
    % Function info text - centered
    uilabel(piece_fig, 'Text', 'Enter any valid MATLAB expression using "t" as the variable.', ...
        'Position', [150, 170, 500, 22]);
    uilabel(piece_fig, 'Text', 'Examples: sin(t), t^2, sqrt(t), pi*cos(t), exp(-t)', ...
        'Position', [150, 150, 500, 22]);
    uilabel(piece_fig, 'Text', 'Note: For section boundaries, make sure end time of section i equals start time of section i+1', ...
        'Position', [150, 130, 500, 22]);
    
    % Buttons centered at bottom
    buttonWidth = 100;
    buttonHeight = 30;
    panelWidth = piece_fig.Position(3);
    buttonY = 60;
    
    % Preview button
    previewButton = uibutton(piece_fig, 'Position', [(panelWidth/2)-150, buttonY, buttonWidth, buttonHeight], ...
        'Text', 'Preview', ...
        'ButtonPushedFcn', @(btn, event) previewPiecewiseFunction(ax, pieceTable));
    
    % Confirm button
    confirmButton = uibutton(piece_fig, 'Position', [(panelWidth/2)-50, buttonY, buttonWidth, buttonHeight], ...
        'Text', 'Confirm', ...
        'ButtonPushedFcn', @(btn, event) confirmPiecewiseFunction(selection_fig, piece_fig, pieceTable));
    
    % Cancel button
    cancelButton = uibutton(piece_fig, 'Position', [(panelWidth/2)+50, buttonY, buttonWidth, buttonHeight], ...
        'Text', 'Cancel', ...
        'ButtonPushedFcn', @(btn, event) cancelInput(piece_fig));
    
    % Handle window close event
    piece_fig.CloseRequestFcn = @(src, event) cancelInput(src);
    
    % Wait for user to confirm or cancel
    uiwait(piece_fig);
end

function previewContinuousFunction(ax, startField, endField, stepField, funcField)
    % Preview the continuous function in the plot area
    try
        % Get time parameters
        start_time = startField.Value;
        end_time = endField.Value;
        time_steps = stepField.Value;
        
        % Validate time parameters
        if end_time <= start_time
            error('End Time must be greater than Start Time.');
        end
        
        if time_steps <= 0
            error('Time Steps must be positive.');
        end
        
        % Get function string
        funcStr = funcField.Value;
        
        % Create a temporary function for preview
        if isempty(funcStr)
            error('Please enter a function expression.');
        end
        
        % Create function handle (using same logic as original get_function_input)
        f = convertStringToFunctionHandle(funcStr);
        
        % Generate time vector and function values
        t = linspace(start_time, end_time, max(round((end_time - start_time) / time_steps) + 1, 2));
        y = arrayfun(f, t);
        
        % Update plot
        cla(ax);
        plot(ax, t, y);
        xlabel(ax, 'Time (s)');
        ylabel(ax, 'r(t)');
        title(ax, sprintf('Function: %s, Time: [%.2f, %.2f]', funcStr, start_time, end_time));
        grid(ax, 'on');
    catch ME
        % Display error message
        uialert(ax.Parent, sprintf('Error: %s', ME.message), 'Error');
    end
end

function previewPiecewiseFunction(ax, pieceTable)
    % Preview the piecewise function in the plot area
    try
        % Get data from table
        data = pieceTable.Data;
        
        % Clear plot
        cla(ax);
        hold(ax, 'on');
        
        % Colors for different sections
        colors = {'b', 'r', 'g', 'm', 'c'};
        
        % Plot each section
        for i = 1:size(data, 1)
            % Get section parameters
            start_time = data{i, 1};
            end_time = data{i, 2};
            time_steps = data{i, 3};
            funcStr = data{i, 4};
            
            % Skip if no function string
            if isempty(funcStr)
                continue;
            end
            
            % Validate time parameters
            if end_time <= start_time
                uialert(ax.Parent, sprintf('Section %d: End Time must be greater than Start Time.', i), 'Error');
                continue;
            end
            
            if time_steps <= 0
                uialert(ax.Parent, sprintf('Section %d: Time Steps must be positive.', i), 'Error');
                continue;
            end
            
            try
                % Create function handle (using same logic as original get_function_input)
                f = convertStringToFunctionHandle(funcStr);
                
                % Generate time vector and function values
                t = linspace(start_time, end_time, max(round((end_time - start_time) / time_steps) + 1, 2));
                y = arrayfun(f, t);
                
                % Plot this section with color coding
                colorIdx = mod(i-1, length(colors)) + 1;
                plot(ax, t, y, colors{colorIdx}, 'LineWidth', 2, 'DisplayName', sprintf('Section %d: %s', i, funcStr));
            catch ME
                % Skip this section if error
                uialert(ax.Parent, sprintf('Error in section %d: %s', i, ME.message), 'Error');
            end
        end
        
        % Finish plot
        legend(ax, 'Location', 'best');
        hold(ax, 'off');
        xlabel(ax, 'Time (s)');
        ylabel(ax, 'r(t)');
        title(ax, 'Piecewise Function Preview');
        grid(ax, 'on');
    catch ME
        % Display error message
        uialert(ax.Parent, sprintf('Error: %s', ME.message), 'Error');
    end
end

function confirmContinuousFunction(selection_fig, cont_fig, startField, endField, stepField, funcField)
    % Confirm the continuous function and store values
    try
        % Get time parameters
        start_time = startField.Value;
        end_time = endField.Value;
        time_steps = stepField.Value;
        
        % Validate time parameters using same logic as get_time_vector
        if end_time <= start_time
            error('End Time must be greater than Start Time.');
        end
        
        if start_time < 0
            error('Start Time must be non-negative.');
        end
        
        if time_steps <= 0
            error('Time Steps must be positive.');
        end
        
        % Check if step size is compatible with time span (validation from get_time_vector)
        numSteps = (end_time - start_time) / time_steps;
        if mod(numSteps, 1) ~= 0
            error('Please enter a suitable size for time steps to divide the time span evenly.');
        end
        
        % Get function string
        funcStr = funcField.Value;
        
        % Create function handle using same format as original get_function_input
        target_value = convertStringToFunctionHandle(funcStr);
        
        % Store values in selection figure's app data
        setappdata(selection_fig, 'target_value', target_value);
        setappdata(selection_fig, 'num_sections', 1);
        setappdata(selection_fig, 'start_time', start_time);
        setappdata(selection_fig, 'end_time', end_time);
        setappdata(selection_fig, 'time_steps', time_steps);
        
        % Close the continuous function figure
        delete(cont_fig);
        
        % Resume execution in main function
        uiresume(selection_fig);
    catch ME
        % Display error message
        uialert(cont_fig, sprintf('Error: %s', ME.message), 'Error');
    end
end

function confirmPiecewiseFunction(selection_fig, piece_fig, pieceTable)
    % Confirm the piecewise function and store values
    try
        % Get data from table
        data = pieceTable.Data;
        num_sections = size(data, 1);
        
        % Extract time vectors
        start_time = zeros(1, num_sections);
        end_time = zeros(1, num_sections);
        time_steps = zeros(1, num_sections);
        
        % Extract and validate time parameters
        for i = 1:num_sections
            start_time(i) = data{i, 1};
            end_time(i) = data{i, 2};
            time_steps(i) = data{i, 3};
            
            % Validate time parameters (same validation as get_time_vector)
            if end_time(i) <= start_time(i)
                error('Section %d: End Time must be greater than Start Time.', i);
            end
            
            if start_time(i) < 0
                error('Section %d: Start Time must be non-negative.', i);
            end
            
            if time_steps(i) <= 0
                error('Section %d: Time Steps must be positive.', i);
            end
            
            % Check if step size is compatible with time span (validation from get_time_vector)
            numSteps = (end_time(i) - start_time(i)) / time_steps(i);
            if mod(numSteps, 1) ~= 0
                error('Section %d: Please enter a suitable size for time steps to divide the time span evenly.', i);
            end
            
            % Check section boundaries
            if i > 1 && abs(start_time(i) - end_time(i-1)) > 1e-10
                error('Section boundaries must match: End time of section %d must equal start time of section %d.', i-1, i);
            end
            
            % Check if function string is empty
            if isempty(data{i, 4})
                error('Function for section %d is empty. Please enter a valid function.', i);
            end
        end
        
        % Create cell array of function handles
        target_value = cell(1, num_sections);
        
        % Convert each function string to a function handle (using same logic as get_function_input)
        for i = 1:num_sections
            funcStr = data{i, 4};
            target_value{i} = convertStringToFunctionHandle(funcStr);
        end
        
        % Store values in selection figure's app data
        setappdata(selection_fig, 'target_value', target_value);
        setappdata(selection_fig, 'num_sections', num_sections);
        setappdata(selection_fig, 'start_time', start_time);
        setappdata(selection_fig, 'end_time', end_time);
        setappdata(selection_fig, 'time_steps', time_steps);
        
        % Close the piecewise function figure
        delete(piece_fig);
        
        % Resume execution in main function
        uiresume(selection_fig);
    catch ME
        % Display error message
        uialert(piece_fig, sprintf('Error: %s', ME.message), 'Error');
    end
end

function f = convertStringToFunctionHandle(funcStr)
    % Convert a string to a function handle using the same logic as get_function_input
    % to maintain compatibility with the original function
    
    % Check if the string is empty
    if isempty(strtrim(funcStr))
        error('Function expression is empty.');
    end
    
    try
        % Check if the string contains 't' (time-dependent)
        if contains(funcStr, 't')
            % Create and test function handle
            f = str2func(['@(t)', funcStr]);
            test_result = f(0); % Test with a dummy value
        else
            % Try to convert to a number or evaluate as a MATLAB expression
            try
                % First try as a number
                eval_value = str2double(funcStr);
                if isnan(eval_value)
                    % Then try as a MATLAB expression
                    eval_value = eval(funcStr);
                end
                % Create function handle for constant value (adding 0*t to make it time-dependent)
                f = str2func(['@(t)', num2str(eval_value), ' + 0*t']);
            catch
                % If all else fails, try directly as a constant function
                f = str2func(['@(t)', funcStr, ' + 0*t']);
            end
        end
    catch ME
        error('Invalid function expression: %s', ME.message);
    end
end

function cancelInput(fig)
    % Handle cancellation
    delete(fig);
end