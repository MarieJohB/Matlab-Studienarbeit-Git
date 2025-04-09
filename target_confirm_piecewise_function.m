function target_confirm_piecewise_function(selection_fig, piece_fig, pieceTable)
% TARGET_CONFIRM_PIECEWISE_FUNCTION - Confirm the piecewise function
%
% Parameters:
%   selection_fig - The main selection figure
%   piece_fig - The piecewise function input figure
%   pieceTable - The table containing section data

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
        % Handle potential string inputs with commas
        if ischar(data{i, 1}) || isstring(data{i, 1})
            start_time(i) = str2double(strrep(data{i, 1}, ',', '.'));
        else
            start_time(i) = data{i, 1};
        end
        
        if ischar(data{i, 2}) || isstring(data{i, 2})
            end_time(i) = str2double(strrep(data{i, 2}, ',', '.'));
        else
            end_time(i) = data{i, 2};
        end
        
        if ischar(data{i, 3}) || isstring(data{i, 3})
            time_steps(i) = str2double(strrep(data{i, 3}, ',', '.'));
        else
            time_steps(i) = data{i, 3};
        end
        
        % Check for NaN values (invalid conversions)
        if isnan(start_time(i)) || isnan(end_time(i)) || isnan(time_steps(i))
            error('Section %d: Please enter valid numeric values for time parameters.', i);
        end
        
        % Validate time parameters
        if end_time(i) <= start_time(i)
            error('Section %d: End Time must be greater than Start Time.', i);
        end
        
        if start_time(i) < 0
            error('Section %d: Start Time must be non-negative.', i);
        end
        
        if time_steps(i) <= 0
            error('Section %d: Time Steps must be positive.', i);
        end
        
        % Check if step size is compatible with time span
        numSteps = (end_time(i) - start_time(i)) / time_steps(i);
        if abs(round(numSteps) - numSteps) > 1e-6
            error('Section %d: Please enter a suitable size for time steps to divide the time span evenly.', i);
        end
        
        % Check section boundaries
        if i > 1 && abs(start_time(i) - end_time(i-1)) > 1e-6
            error('Section boundaries must match: End time of section %d must equal start time of section %d.', i-1, i);
        end
        
        % Check if function string is empty
        if isempty(data{i, 4})
            error('Function for section %d is empty. Please enter a valid function.', i);
        end
    end
    
    % Create cell array of function handles
    target_value = cell(1, num_sections);
    
    % Convert each function string to a function handle
    for i = 1:num_sections
        funcStr = data{i, 4};
        % Replace commas with periods in function string
        funcStr = strrep(funcStr, ',', '.');
        target_value{i} = convert_string_to_function_handle(funcStr);
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
    uialert(piece_fig, ME.message, 'Error');
end
end