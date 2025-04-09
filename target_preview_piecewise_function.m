function target_preview_piecewise_function(ax, pieceTable)
% TARGET_PREVIEW_PIECEWISE_FUNCTION - Preview the piecewise function in the plot area
%
% Parameters:
%   ax - The axes to plot in
%   pieceTable - The table containing section data

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
        % Get section parameters and handle potential comma as decimal separator
        if ischar(data{i, 1}) || isstring(data{i, 1})
            start_time = str2double(strrep(data{i, 1}, ',', '.'));
        else
            start_time = data{i, 1};
        end
        
        if ischar(data{i, 2}) || isstring(data{i, 2})
            end_time = str2double(strrep(data{i, 2}, ',', '.'));
        else
            end_time = data{i, 2};
        end
        
        if ischar(data{i, 3}) || isstring(data{i, 3})
            time_steps = str2double(strrep(data{i, 3}, ',', '.'));
        else
            time_steps = data{i, 3};
        end
        
        funcStr = data{i, 4};
        
        % Skip if no function string
        if isempty(funcStr)
            continue;
        end
        
        % Validate time parameters
        if isnan(start_time) || isnan(end_time) || isnan(time_steps)
            uialert(ax.Parent, sprintf('Section %d: Invalid numeric values.', i), 'Error');
            continue;
        end
        
        if end_time <= start_time
            uialert(ax.Parent, sprintf('Section %d: End Time must be greater than Start Time.', i), 'Error');
            continue;
        end
        
        if time_steps <= 0
            uialert(ax.Parent, sprintf('Section %d: Time Steps must be positive.', i), 'Error');
            continue;
        end
        
        try
            % Create function handle
            % Replace comma with period in function string
            funcStr = strrep(funcStr, ',', '.');
            f = convert_string_to_function_handle(funcStr);
            
            % Generate time vector and function values
            t = linspace(start_time, end_time, max(round((end_time - start_time) / time_steps) + 1, 2));
            y = arrayfun(f, t);
            
            % Plot this section with color coding
            colorIdx = mod(i-1, length(colors)) + 1;
            plot(ax, t, y, colors{colorIdx}, 'DisplayName', sprintf('Section %d: %s', i, funcStr));
        catch ME
            % Skip this section if error
            uialert(ax.Parent, sprintf('Error in section %d: %s', i, ME.message), 'Error');
        end
    end
    
    % Finish plot
    legend(ax, 'Location', 'best');
    hold(ax, 'off');
    xlabel(ax, 'Time (s)');
    % No y-label as requested
    title(ax, 'Piecewise Function Preview');
    grid(ax, 'on');
catch ME
    % Display error message
    cla(ax);
    text(ax, 0.5, 0.5, ['Error: ', ME.message], ...
         'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'red');
end
end