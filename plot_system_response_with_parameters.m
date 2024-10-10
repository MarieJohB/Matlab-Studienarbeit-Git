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
    yline(a, 'k--', 'DisplayName', 'Stationary Value (a)');
    yline(a + b, 'g--', 'DisplayName', 'Overshoot (b)');
    yline(a - c, 'm--', 'DisplayName', 'Decay Ratio (c)');
    xline(tr, 'c--', 'DisplayName', 'Rise Time (tr)');
    xline(ts, 'y--', 'DisplayName', 'Settling Time (ts)');
    title('System Response with Key Parameters');
    xlabel('Time (seconds)');
    ylabel('Response');
    legend;

    % Annotate the plot with calculated values
    annotation('textbox', [0.15, 0.7, 0.2, 0.1], 'String', ['Stationary value (a): ', num2str(a)], 'EdgeColor', 'none', 'FontSize', 12, 'Color', 'black');
    annotation('textbox', [0.15, 0.65, 0.2, 0.1], 'String', ['Overshoot (b): ', num2str(b)], 'EdgeColor', 'none', 'FontSize', 12, 'Color', 'black');
    annotation('textbox', [0.15, 0.6, 0.2, 0.1], 'String', ['Decay ratio (c): ', num2str(c)], 'EdgeColor', 'none', 'FontSize', 12, 'Color', 'black');
    annotation('textbox', [0.15, 0.55, 0.2, 0.1], 'String', ['Rise time (tr): ', num2str(tr)], 'EdgeColor', 'none', 'FontSize', 12, 'Color', 'black');
    annotation('textbox', [0.15, 0.5, 0.2, 0.1], 'String', ['Settling time (ts): ', num2str(ts)], 'EdgeColor', 'none', 'FontSize', 12, 'Color', 'black');

    grid on;

    % Display key parameters in the command window
    disp(['Stationary value (a): ', num2str(a)]);
    disp(['Overshoot (b): ', num2str(b)]);
    disp(['Decay ratio (c): ', num2str(c)]);
    disp(['Rise time (tr): ', num2str(tr)]);
    disp(['Settling time (ts): ', num2str(ts)]);
end
