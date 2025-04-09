function target_confirm_continuous_function(selection_fig, cont_fig, startField, endField, stepField, funcField)
try
    % Get time parameters
    start_time = startField.Value;
    end_time = endField.Value;
    time_steps = stepField.Value;
    
    % Validate time parameters
    if end_time <= start_time
        error('End Time must be greater than Start Time.');
    end
    
    if start_time < 0
        error('Start Time must be non-negative.');
    end
    
    if time_steps <= 0
        error('Time Steps must be positive.');
    end
    
    % Check if step size is compatible with time span
    numSteps = (end_time - start_time) / time_steps;
    if mod(numSteps, 1) ~= 0
        error('Please enter a suitable size for time steps to divide the time span evenly.');
    end
    
    % Get function string and replace commas with periods
    funcStr = strrep(funcField.Value, ',', '.');
    
    % Create function handle
    target_value = convert_string_to_function_handle(funcStr);
    
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
    uialert(cont_fig, ME.message, 'Error');
end
end