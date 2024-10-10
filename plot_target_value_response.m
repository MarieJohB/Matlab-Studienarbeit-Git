function plot_target_value_response(G, target_value)
    % Display the transfer function
    disp('Transfer Function:');
    disp(G);

    % Define the time vector
    t = linspace(0, 10, 1000); % 0 to 10 seconds with 1000 points

    % Calculate the system response based on the type of target value
    if isa(target_value, 'function_handle')
        target = arrayfun(target_value, t); % Evaluate the function handle over time
    else
        target = target_value * ones(size(t)); % Constant target value
    end

    % Define the unity feedback system
    closed_loop_system = feedback(G, 1);

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