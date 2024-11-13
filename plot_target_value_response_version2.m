function plot_target_value_response_version2(L, target_value, start_time, end_time, time_steps)
    % first: check whether the target value was defined piecewise or continuous
    num_sections = length(target_value);

    if num_sections == 1 % target value is continuous
        % Define the time vector
        num_steps = (end_time - start_time) / time_steps; % calculate number of steps
        t = linspace(start_time, end_time, num_steps + 1)'; % create time vector with specified values, transpose to column vector
        u = target_value;

        % Calculate the system response based on the type of target value
        if isa(u, 'function_handle')
            target = arrayfun(u, t); % Evaluate the function handle over time
        else
            target = u * ones(size(t)); % Constant target value
        end

    elseif num_sections > 1 % target value is defined piecewise
        % Initialize the time vector and target array
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

            if isa(target_value{i}, 'function_handle')
                % Evaluate the function handle over the current time section
                func_target = arrayfun(target_value{i}, t_section);
                targets = [targets; func_target]; % Append as new rows
            else
                % Constant target value over the current time section
                const_target = repmat(target_value{i}, num_steps + 1, 1);
                targets = [targets; const_target]; % Append as new rows
            end
        end
        
        target = targets; % Ensure correct dimensions

    else
        disp('Something went wrong. Number of sections is smaller than 1.');
        return; % End the function if num_sections is smaller than 1
    end

    % Debug-Ausgaben
    disp('L:');
    disp(L);
    disp('Target:');
    disp(size(target)); % Output the size of target
    disp('Time vector t:');
    disp(size(t)); % Output the size of t
    disp('Time vector t values:');
    disp(t'); % Display the actual time vector values for verification

    % Ensure the time vector is monotonically increasing and evenly spaced
    if any(diff(t) <= 0)
        error('The time vector t must be monotonically increasing and evenly spaced.');
    end

    % Define the unity feedback system
    closed_loop_system = feedback(L, 1); % Ensure both arguments are passed

    % Calculate the system response
    [y, t_out] = lsim(closed_loop_system, target, t);

    % Plot the target value and system response
    figure;
    plot(t, target, 'r--', 'DisplayName', 'Target Value');
    hold on;
    plot(t_out, y, 'b', 'DisplayName', 'System Response');
    title('System Response to Target Value');
    xlabel('Time (seconds)');
    ylabel('Amplitude');
    legend;
    grid on;

    % Display the closed-loop transfer function
    disp('Closed-Loop Transfer Function:');
    disp(closed_loop_system);
end
