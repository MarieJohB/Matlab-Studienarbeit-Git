function plot_target_value_response_version2(G, K, target_value)
    % first: check whether the target value was defined piecewise or continuous
    num_sections = length(target_value);

    % call function to get time sections
    [start_time, end_time, time_steps] = get_time_vector(num_sections);
  

    if num_sections == 1 % target value is continuous
        % Define the time vector
        [t] = create_linear_time_vector(1);
        u = target_value;

       
        if isa(u, 'function_handle')
            target = arrayfun(u, t); % Evaluate the function handle over time
        else
            target = u * ones(size(t)); % Constant target value
        end

    elseif num_sections > 1 % target value is defined piecewise
        
        [t] = create_linear_time_vector(num_sections);
        [values_vector] = generate_Function_Values_Version2(target_value);
        target = values_vector;

    else
        disp('Something went wrong. Number of sections is smaller than 1.');
        return; % End the function if num_sections is smaller than 1
    end

    % Debug-Ausgaben
    disp('Target:');
    disp(size(target)); % Output the size of target
    disp('Time vector t:');
    disp(size(t)); % Output the size of t
    
    disp(target); % test purposes

    % Ensure the time vector is monotonically increasing and evenly spaced
    if any(diff(t) <= 0)
        error('The time vector t must be monotonically increasing and evenly spaced.');
    end

    
    [T, S, L, GS, KS] = transferfunctions(G, K);



    % Calculate the system response
    [y, t_out] = lsim(T, target, t);

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


end
