function disturbance_preview_continuous_function(ax, startField, endField, stepField, funcField, disturbance_label)
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
    
    % Replace commas with periods for numerical consistency
    funcStr = strrep(funcStr, ',', '.');
    
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
    title(ax, ['Function: ', funcStr]);
    grid(ax, 'on');
catch ME
    % Display error message
    cla(ax);
    text(ax, 0.5, 0.5, ['Error: ', ME.message], ...
         'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'red');
end
end