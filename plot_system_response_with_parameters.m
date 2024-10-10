function plot_system_response_with_parameters(L, target_value, d1, d2)
    % Display the transfer function
    disp('Transfer Function:');
    disp(L);

    % Get the timeframe for the plot from the user
    prompt = {'Enter the start time (seconds):', 'Enter the end time (seconds):'};
    dlgtitle = 'Timeframe Input';
    dims = [1 50];
    answer = inputdlg(prompt, dlgtitle, dims);
    
    % Validate input
    try
        start_time = str2double(answer{1});
        end_time = str2double(answer{2});
        if isnan(start_time) || isnan(end_time) || start_time < 0 || end_time <= start_time
            error('Invalid input. Please enter valid time values.');
        end
    catch
        uiwait(msgbox('Invalid input. Please enter the values correctly.', 'Error', 'error'));
        return;
    end

    % Define the time vector
    t = linspace(start_time, end_time, 1000); % Start time to end time with 1000 points

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
    c = a - min(y_total); % Decay ratio
    tr = t(find(y_total >= a * 0.9, 1)); % Rise time (90% of a)
    ts = t(find(abs(y_total - a) <= 0.02 * a, 1, 'last')); % Settling time (2% band)

    % Plot the system response with key parameters
    figure;
    plot(t_out, y_total, 'b', 'DisplayName', 'System Response');
    hold on;
    plot(t_out, target, 'r--', 'DisplayName', ['Target Value: ', num2str(target_value), ' units']);
    yline(a, 'k--', 'DisplayName', ['Stationary Value (a): ', num2str(a), ' units']);
    yline(a + b, 'g--', 'DisplayName', ['Overshoot (b): ', num2str(b), ' units']);
    yline(a - c, 'm--', 'DisplayName', ['Decay Ratio (c): ', num2str(c), ' units']);
    xline(tr, 'c--', 'DisplayName', ['Rise Time (tr): ', num2str(tr), ' seconds']);
    xline(ts, 'y--', 'DisplayName', ['Settling Time (ts): ', num2str(ts), ' seconds']);
    title('System Response with Key Parameters');
    xlabel('Time (seconds)');
    ylabel('Response (units)');
    legend('Location', 'best');
    grid on;

    % Display key parameters in the command window
    disp(['Stationary value (a): ', num2str(a), ' units']);
    disp(['Overshoot (b): ', num2str(b), ' units']);
    disp(['Decay ratio (c): ', num2str(c), ' units']);
    disp(['Rise time (tr): ', num2str(tr), ' seconds']);
    disp(['Settling time (ts): ', num2str(ts), ' seconds']);
end

