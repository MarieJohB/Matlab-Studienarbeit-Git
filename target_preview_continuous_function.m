function target_preview_continuous_function(ax, startField, endField, stepField, funcField)
try
    % Get time parameters and convert to numbers
    start_time = str2double(strrep(startField.Value, ',', '.'));
    end_time = str2double(strrep(endField.Value, ',', '.'));
    time_steps = str2double(strrep(stepField.Value, ',', '.'));
    
    % Validate time parameters
    if isnan(start_time) || isnan(end_time) || isnan(time_steps)
        error('Please enter valid numeric values for time parameters.');
    end
    
    if end_time <= start_time
        error('End Time must be greater than Start Time.');
    end
    
    if time_steps <= 0
        error('Time Steps must be positive.');
    end
    
    % Get function string
    funcStr = funcField.Value;
    
    % If function is empty, show empty plot
    if isempty(funcStr)
        cla(ax);
        title(ax, 'Function Preview');
        return;
    end
    
    % Create a temporary function for preview
    f = convert_string_to_function_handle(funcStr);
    
    % Generate time vector and function values
    t = linspace(start_time, end_time, max(round((end_time - start_time) / time_steps) + 1, 2));
    y = arrayfun(f, t);
    
    % Update plot
    cla(ax);
    plot(ax, t, y);
    xlabel(ax, 'Time (s)');
    % No y-label as requested
    title(ax, ['Function: ', funcStr]);
    grid(ax, 'on');
catch ME
    % Display error message
    cla(ax);
    text(ax, 0.5, 0.5, ['Error: ', ME.message], ...
         'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'red');
end
end