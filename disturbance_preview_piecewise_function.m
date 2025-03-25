function disturbance_preview_piecewise_function(ax, pieceTable, disturbance_label)
% DISTURBANCE_PREVIEW_PIECEWISE_FUNCTION - Preview the piecewise function in the plot area
%
% Parameters:
%   ax - The axes to plot in
%   pieceTable - The table containing section data
%   disturbance_label - The label for the disturbance (e.g. 'd_1')

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
            % Create function handle
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