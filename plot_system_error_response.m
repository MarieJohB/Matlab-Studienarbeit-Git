function plot_system_error_response(L, target_value, d1, d2)
    % Display the transfer function
    disp('Transfer Function:');
    disp(L);

    % Define the time vector
    t = linspace(0, 10, 1000); % 0 to 10 seconds with 1000 points

    % Calculate the target value
    if isa(target_value, 'function_handle')
        target = arrayfun(target_value, t); % Evaluate the function handle over time
    else
        target = target_value * ones(size(t)); % Constant target value
    end

    % Define the unity feedback system
    closed_loop_system = feedback(L, 1);

    % Calculate the system response to the target value
    [y_target, t_out] = lsim(closed_loop_system, target, t);

    % Calculate the system response to disturbances
    H_d1 = feedback(1, L); % Transfer function for disturbance D1
    H_d2 = L * feedback(1, L); % Transfer function for disturbance D2

    [y_d1, ~] = lsim(H_d1, d1, t);
    [y_d2, ~] = lsim(H_d2, d2, t);

    % Total system response
    y_total = y_target + y_d1 + y_d2;

    % Calculate the error
    error = target - y_total;

    % Plot the error response
    figure;
    plot(t_out, error, 'r', 'DisplayName', 'System Error');
    title('System Error Response');
    xlabel('Time (seconds)');
    ylabel('Error');
    legend('System Error'); % Ensure legend displays correctly
    grid on;

    % Display the closed-loop transfer function
    disp('Closed-Loop Transfer Function:');
    disp(closed_loop_system);
end
