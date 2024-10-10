function get_step_response(T)
    % Display the transfer function
    disp('Transfer Function:');
    disp(T);

    % Define the unity feedback system
    closed_loop_system = feedback(T, 1);

    % Calculate and plot the step response
    figure;
    step(closed_loop_system);
    title('Closed-Loop Step Response');
    xlabel('Time (seconds)');
    ylabel('Amplitude');
    grid on;

    % Display the closed-loop transfer function
    disp('Closed-Loop Transfer Function:');
    disp(closed_loop_system);
end