function [values_vector] = generate_Function_Values_Version2(function_vector)
    % Creates values of a time dependent function 
    % Input can be a function that was defined continuous or piecewise
    % Output: vector with values that can be plotted with corresponding time
    % vector t 
    % Input target_values -> Output: values of r(t)
    % or input is disturbance -> output: values of d(t)

    % First: check whether the function was defined piecewise or continuous
    num_sections = length(function_vector);

    % Call function to get the time vector
    [t] = create_linear_time_vector(num_sections);

    % Initialize the values vector
    values_vector = zeros(size(t));

    % Call function to get time sections
    [start_time, end_time, time_steps] = get_time_vector(num_sections);

    % Validate the time sections
   % if length(start_time) ~= num_sections || length(end_time) ~= num_sections || length(time_steps) ~= num_sections
    %    error('The lengths of start_time, end_time, and time_steps must match num_sections.');
    %end

    if num_sections == 1 % Function is continuous
        u = function_vector;

        if isa(u, 'function_handle')
            values_vector = arrayfun(u, t); % Evaluate the function handle over time
        else
            values_vector = u * ones(size(t)); % Constant function
        end

    elseif num_sections > 1 % Function is defined piecewise
        % Initialize index for target values
        idx = 1;

        for i = 1:num_sections
            % Determine the time section
            t_section_start = find(t >= start_time(i), 1, 'first');
            t_section_end = find(t <= end_time(i), 1, 'last');
            t_section = t(t_section_start:t_section_end);

            % Avoid duplicate values at the transition points
            if i > 1 && t_section(1) == t(idx-1)
                t_section = t_section(2:end); % Remove the first element to avoid duplicate
            end

            if isa(function_vector{i}, 'function_handle')
                % Evaluate the function handle over the current time section
                func_target = arrayfun(function_vector{i}, t_section);
                values_vector(idx:idx+length(func_target)-1) = func_target;
            else
                % Constant function over the current time section
                const_target = repmat(function_vector{i}, length(t_section), 1);
                values_vector(idx:idx+length(const_target)-1) = const_target;
            end

            % Update the index
            idx = idx + length(t_section);
        end

    else
        disp('Something went wrong. Number of sections is smaller than 1.');
        return; % End the function if num_sections is smaller than 1
    end

    % Debugging: Check and display dimensions of t and values_vector
    disp('Debugging Info:');
    disp(['Length of t: ', num2str(length(t))]);
    disp(['Length of values_vector: ', num2str(length(values_vector))]);

    if length(t) ~= length(values_vector)
        disp('Error: The lengths of t and values_vector do not match!');
    else
        disp('Success: The lengths of t and values_vector match.');
    end
end
