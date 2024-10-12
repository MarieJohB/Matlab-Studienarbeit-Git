function K = get_user_controller()
    % Controller selection dialog
    controller_type = menu('Choose a controller type:', ...
        'Proportional (P)', 'Proportional-Integral (PI)', 'Proportional-Derivative (PD)', ...
        'Proportional-Integral-Derivative (PID)', 'PT1', 'PIT1', 'I2', 'PIDT1', 'Custom');
    
    % Initialize coefficients
    Kp = 0; Ki = 0; Kd = 0; t_controller = 0;
    
    % Define prompt and dialog title based on controller type
    switch controller_type
        case 1 % P Controller
            prompt = {'Enter the proportional gain Kp:'}; dlgtitle = 'P Controller Input';
        case 2 % PI Controller
            prompt = {'Enter the proportional gain Kp:', 'Enter the integral gain Ki:'}; dlgtitle = 'PI Controller Input';
        case 3 % PD Controller
            prompt = {'Enter the proportional gain Kp:', 'Enter the derivative gain Kd:'}; dlgtitle = 'PD Controller Input';
        case 4 % PID Controller
            prompt = {'Enter the proportional gain Kp:', 'Enter the integral gain Ki:', 'Enter the derivative gain Kd:'}; dlgtitle = 'PID Controller Input';
        case 5 % PT1 Controller
            prompt = {'Enter the proportional gain Kp:', 'Enter the time constant T:'}; dlgtitle = 'PT1 Controller Input';
        case 6 % PIT1 Controller
            prompt = {'Enter the proportional gain Kp:', 'Enter the integral gain Ki:', 'Enter the time constant T:'}; dlgtitle = 'PIT1 Controller Input';
        case 7 % I2 Controller
            prompt = {'Enter the integral gain Ki:'}; dlgtitle = 'I2 Controller Input';
        case 8 % PIDT1 Controller
            prompt = {'Enter the proportional gain Kp:', 'Enter the integral gain Ki:', 'Enter the derivative gain Kd:', 'Enter the time constant T:'}; dlgtitle = 'PIDT1 Controller Input';
        case 9 % Custom Controller
            prompt = {'Enter the numerator coefficients as a vector [b0 b1 b2 ...]:', ...
                      'Enter the denominator coefficients as a vector [a0 a1 a2 ...]:'}; dlgtitle = 'Custom Controller Input';
        otherwise
            disp('No controller selected.'); K = []; return;
    end
    
    % Input dimensions
    dims = [1 50];
    
    % Initialize K to be empty
    K = [];
    
    % Loop until valid input is provided or user cancels the dialog
    while isempty(K)
        % Display input dialog
        answer = inputdlg(prompt, dlgtitle, dims);

        % Check if user cancelled the dialog
        if isempty(answer)
            disp('Operation cancelled by user.');
            return;
        end

        % Validate and assign inputs based on controller type
        try
            switch controller_type
                case 1 % P Controller
                    Kp_str = strrep(answer{1}, ',', '.'); % checking for "," and replacing with "."
                    Kp = str2double(Kp_str);
                    K = pid(Kp, 0, 0);
                case 2 % PI Controller
                    Kp_str = strrep(answer{1}, ',', '.'); % checking for "," and replacing with "."
                    Kp = str2double(Kp_str);
                    Ki_str = strrep(answer{2}, ',', '.'); % checking for "," and replacing with "."
                    Ki = str2double(Ki_str);
                    K = pid(Kp, Ki, 0);
                case 3 % PD Controller
                    Kp_str = strrep(answer{1}, ',', '.'); % checking for "," and replacing with "."
                    Kp = str2double(Kp_str);
                    Kd_str = strrep(answer{2}, ',', '.'); % checking for "," and replacing with "."
                    Kd = str2double(Kd_str);
                    K = pid(Kp, 0, Kd);
                case 4 % PID Controller
                    Kp_str = strrep(answer{1}, ',', '.'); % checking for "," and replacing with "."
                    Kp = str2double(Kp_str);
                    Ki_str = strrep(answer{2}, ',', '.'); % checking for "," and replacing with "."
                    Ki = str2double(Ki_str);
                    Kd_str = strrep(answer{3}, ',', '.'); % checking for "," and replacing with "."
                    Kd = str2double(Kd_str);
                    K = pid(Kp, Ki, Kd);
                case 5 % PT1 Controller
                    Kp_str = strrep(answer{1}, ',', '.'); % checking for "," and replacing with "."
                    Kp = str2double(Kp_str);
                    t_controller_str = strrep(answer{2}, ',', '.'); % checking for "," and replacing with "."
                    t_controller = str2double(t_controller_str);
                    K = tf([Kp], [t_controller 1]);
                case 6 % PIT1 Controller
                    Kp_str = strrep(answer{1}, ',', '.'); % checking for "," and replacing with "."
                    Kp = str2double(Kp_str);
                    Ki_str = strrep(answer{2}, ',', '.'); % checking for "," and replacing with "."
                    Ki = str2double(Ki_str);
                    t_controller_str = strrep(answer{3}, ',', '.'); % checking for "," and replacing with "."
                    t_controller = str2double(t_controller_str);
                    K = tf([Kp*Ki Kp], [t_controller 1 0]);
                case 7 % I2 Controller
                    Ki_str = strrep(answer{1}, ',', '.'); % checking for "," and replacing with "."
                    Ki = str2double(Ki_str);
                    K = tf([Ki], [1 0 0]);
                case 8 % PIDT1 Controller
                    Kp_str = strrep(answer{1}, ',', '.'); % checking for "," and replacing with "."
                    Kp = str2double(Kp_str);
                    Ki_str = strrep(answer{2}, ',', '.'); % checking for "," and replacing with "."
                    Ki = str2double(Ki_str);
                    Kd_str = strrep(answer{3}, ',', '.'); % checking for "," and replacing with "."
                    Kd = str2double(Kd_str);
                    t_controller_str = strrep(answer{4}, ',', '.'); % checking for "," and replacing with "."
                    t_controller = str2double(t_controller_str);
                    K = tf([Kd Kp Ki], [t_controller 1 0]);
                case 9 % Custom Controller
                    num_str = strrep(answer{1}, ',', '.'); % checking for "," and replacing with "."
                    num = str2num(num_str); %#ok<ST2NM>
                    den_str = strrep(answer{2}, ',', '.'); % checking for "," and replacing with "."
                    den = str2num(den_str); %#ok<ST2NM>
                    K = tf(num, den);
            end
        catch
            uiwait(msgbox('Invalid input. Please enter the gains and time constant correctly.', 'Error', 'error'));
            K = [];
        end
    end

    % Display the controller transfer function
    if ~isempty(K)
        disp('The entered controller is:');
        disp(K);
    end
end
