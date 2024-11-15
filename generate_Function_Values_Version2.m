function [values_vector] = generate_Function_Values_Version2(function_vector)
% creates values of a time dependent function 
% input can be a function that was defined continous or piecewise
% output: vector with values that can be plotted with corresponding time
% vector t 
% input traget_values -> Output: values of r(t)
% or input is disturbance -> output: values of d(t)



    % first: check whether the function was defined piecewise or continuous
    num_sections = length(function_vector);

    % call function to get time sections
    [start_time, end_time, time_steps] = get_time_vector(num_sections);
  

    if num_sections == 1 % function is continuous
        % Define the time vector
        [t] = create_linear_time_vector(1);
        u = function_vector;

       
        if isa(u, 'function_handle')
            values_vector = arrayfun(u, t); % Evaluate the function handle over time
        else
            values_vector = u * ones(size(t)); % Constant function
        end

    elseif num_sections > 1 % function is defined piecewise
        % Initialize the time vector and function array
        t = [];
        targets = [];
        
        for i = 1:num_sections
            % useing min(time_steps) for all because when using lsim: time has to
            % be spaced venly
            num_steps = round((end_time(i) - start_time(i)) / min(time_steps)); % calculate integer number of steps
            t_section = linspace(start_time(i), end_time(i), num_steps + 1)'; % create equal intervals for each section with time_steps(i)
            
            % Avoid duplicate values at the transition points
            % neccessary because end_time(1) = start_time(2) and so on 
            if i > 1
                t_section = t_section(2:end); % remove the first element to avoid duplicate
            end

            t = [t; t_section]; % concatenate as column vector

            if isa(function_vector{i}, 'function_handle')
                % Evaluate the function handle over the current time section
                func_target = arrayfun(function_vector{i}, t_section);
                targets = [targets; func_target]; % Append as new rows
            else
                % Constant function over the current time section
                const_target = repmat(function_vector{i}, num_steps + 1, 1);
                targets = [targets; const_target]; % Append as new rows
            end
        end
        
        values_vector = targets; % Ensure correct dimensions

    else
        disp('Something went wrong. Number of sections is smaller than 1.');
        return; % End the function if num_sections is smaller than 1
    end


end

