function plot_system_response_with_parameters(L, target_value, d1, d2)
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

    % Calculate key parameters
    a = y_total(end); % Stationary value
    b = max(y_total) - a; % Overshoot
    c = b - min(y_total); % Decay ratio
    tr = t(find(y_total >= a * 0.9, 1)); % Rise time (90% of a)
    ts = t(find(abs(y_total - a) <= 0.02 * a, 1, 'last')); % Settling time (2% band)

    % Plot the system response with key parameters
    figure;
    plot(t_out, y_total, 'b', 'DisplayName', 'System Response');
    hold on;
    plot(t_out, target, 'r--', 'DisplayName', 'Target Value');
    plot([0 t_out(end)], [a a], 'k--', 'DisplayName', 'Stationary Value (a)');
    plot([0 t_out(end)], [a + b a + b], 'g--', 'DisplayName', 'Overshoot (b)');
    plot([0 t_out(end)], [a - c a - c], 'm--', 'DisplayName', 'Decay Ratio (c)');
    plot([tr tr], ylim, 'c--', 'DisplayName', 'Rise Time (tr)');
    plot([ts ts], ylim, 'y--', 'DisplayName', 'Settling Time (ts)');
    title('System Response with Key Parameters');
    xlabel('Time (seconds)');
    ylabel('Response');
    legend;
    grid on;

    % Display key parameters
    disp(['Stationary value (a): ', num2str(a)]);
    disp(['Overshoot (b): ', num2str(b)]);
    disp(['Decay ratio (c): ', num2str(c)]);
    disp(['Rise time (tr): ', num2str(tr)]);
    disp(['Settling time (ts): ', num2str(ts)]);
end