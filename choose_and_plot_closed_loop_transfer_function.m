function choose_and_plot_closed_loop_transfer_function(L)
    % Transfer function selection dialog
    TF_Type = menu('Choose a closed-loop transfer function to plot:', ...
                   'Y/R (Input to Output) ', ...
                   'Y/H (Disturbance to Output)', ...
                   'Y/D1 (Disturbance 1 to Output) ', ...
                   'Y/D2 (Disturbance 2 to Output)', ...
                   'U/R (Input to Control Effort)', ...
                   'U/H (Disturbance to Control Effort)', ...
                   'E/R (Input to Error)', ...
                   'E/H (Disturbance to Error)', ...
                   'E/D1 (Disturbance 1 to Error) ', ...
                   'E/D2 (Disturbance 2 to Error) ', ...
                   'U/D1 (Disturbance 1 to Control Effort) ', ...
                   'U/D2 (Disturbance 2 to Control Effort) ');

    % Initialize closed-loop transfer function
    H = [];

    % Define transfer function based on user choice
    switch TF_Type
        case 1 % Y/R (Input to Output) - G_{yr}(s)
            H = feedback(L, 1);
        case 2 % Y/H (Disturbance to Output) - G_{yh}(s)
            H = -feedback(L, 1);
        case 3 % Y/D1 (Disturbance 1 to Output) - G_{yd1}(s)
            H = feedback(1, L);
        case 4 % Y/D2 (Disturbance 2 to Output) - G_{yd2}(s)
            H = L * feedback(1, L);
        case 5 % U/R (Input to Control Effort) - G_{ur}(s)
            H = L * feedback(1, L);
        case 6 % U/H (Disturbance to Control Effort) - G_{uh}(s)
            H = -L * feedback(1, L);
        case 7 % E/R (Input to Error) - G_{er}(s)
            H = feedback(1, L);
        case 8 % E/H (Disturbance to Error) - G_{eh}(s)
            H = -feedback(1, L);
        case 9 % E/D1 (Disturbance 1 to Error) - G_{ed1}(s)
            H = -feedback(1, L);
        case 10 % E/D2 (Disturbance 2 to Error) - G_{ed2}(s)
            H = -L * feedback(1, L);
        case 11 % U/D1 (Disturbance 1 to Control Effort) - G_{ud1}(s)
            H = -L * feedback(1, L);
        case 12 % U/D2 (Disturbance 2 to Control Effort) - G_{ud2}(s)
            H = -feedback(L, 1);
        otherwise
            disp('No valid transfer function selected.');
            return;
    end

    % Plot step response if a valid transfer function is chosen
    if ~isempty(H)
        figure;
        step(H);
        title('Step Response');
        xlabel('Time (seconds)');
        ylabel('Amplitude');
        grid on;
        disp('Selected Transfer Function:');
        disp(H);
    else
        disp('No valid transfer function to plot.');
    end
end
