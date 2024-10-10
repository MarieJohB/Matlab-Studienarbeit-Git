function K = get_user_controller()
    % Controller selection dialog
    controller_type = menu('Choose a controller type:', ...
        'Proportional (P)', 'Proportional-Integral (PI)', 'Proportional-Derivative (PD)', ...
        'Proportional-Integral-Derivative (PID)', 'PT1', 'PIT1', 'I2', 'PIDT1', 'Custom');

    % Initialize coefficients
    Kp = 0;
    Ki = 0;
    Kd = 0;
    t_controller = 0;

    switch controller_type
        case 1 % P Controller
            prompt = {'Enter the proportional gain Kp:'};
            dlgtitle = 'P Controller Input';
        case 2 % PI Controller
            prompt = {'Enter the proportional gain Kp:', 'Enter the integral gain Ki:'};
            dlgtitle = 'PI Controller Input';
        case 3 % PD Controller
            prompt = {'Enter the proportional gain Kp:', 'Enter the derivative gain Kd:'};
            dlgtitle = 'PD Controller Input';
        case 4 % PID Controller
            prompt = {'Enter the proportional gain Kp:', 'Enter the integral gain Ki:', 'Enter the derivative gain Kd:'};
            dlgtitle = 'PID Controller Input';
        case 5 % PT1 Controller
            prompt = {'Enter the proportional gain Kp:', 'Enter the time constant T:'};
            dlgtitle = 'PT1 Controller Input';
        case 6 % PIT1 Controller
            prompt = {'Enter the proportional gain Kp:', 'Enter the integral gain Ki:', 'Enter the time constant T:'};
            dlgtitle = 'PIT1 Controller Input';
        case 7 % I2 Controller
            prompt = {'Enter the integral gain Ki:'};
            dlgtitle = 'I2 Controller Input';
        case 8 % PIDT1 Controller
            prompt = {'Enter the proportional gain Kp:', 'Enter the integral gain Ki:', 'Enter the derivative gain Kd:', 'Enter the time constant T:'};
            dlgtitle = 'PIDT1 Controller Input';
        case 9 % Custom Controller
            prompt = {'Enter the numerator coefficients as a vector [b0 b1 b2 ...]:', ...
                      'Enter the denominator coefficients as a vector [a0 a1 a2 ...]:'};
            dlgtitle = 'Custom Controller Input';
        otherwise
            disp('No controller selected.');
            K = [];
            return;
    end

    % Input dimensions
    dims = [1 50];

    % Display input dialog
    answer = inputdlg(prompt, dlgtitle, dims);

    % Validate and assign inputs based on controller type
    try
        switch controller_type
            case 1 % P Controller
                Kp = str2double(answer{1});
                K = pid(Kp, 0, 0);
            case 2 % PI Controller
                Kp = str2double(answer{1});
                Ki = str2double(answer{2});
                K = pid(Kp, Ki, 0);
            case 3 % PD Controller
                Kp = str2double(answer{1});
                Kd = str2double(answer{2});
                K = pid(Kp, 0, Kd);
            case 4 % PID Controller
                Kp = str2double(answer{1});
                Ki = str2double(answer{2});
                Kd = str2double(answer{3});
                K = pid(Kp, Ki, Kd);
            case 5 % PT1 Controller
                Kp = str2double(answer{1});
                t_controller = str2double(answer{2});
                K = tf([Kp], [t_controller 1]);
            case 6 % PIT1 Controller
                Kp = str2double(answer{1});
                Ki = str2double(answer{2});
                t_controller = str2double(answer{3});
                K = tf([Kp*Ki Kp], [t_controller 1 0]);
            case 7 % I2 Controller
                Ki = str2double(answer{1});
                K = tf([Ki], [1 0 0]);
            case 8 % PIDT1 Controller
                Kp = str2double(answer{1});
                Ki = str2double(answer{2});
                Kd = str2double(answer{3});
                t_controller = str2double(answer{4});
                K = tf([Kd Kp Ki], [t_controller 1 0]);
            case 9 % Custom Controller
                num = str2num(answer{1}); %#ok<ST2NM>
                den = str2num(answer{2}); %#ok<ST2NM>
                K = tf(num, den);
        end
    catch
        uiwait(msgbox('Invalid input. Please enter the gains and time constant correctly.', 'Error','error'));
        K = [];
    end

    % Display the controller transfer function
    if ~isempty(K)
        disp('The entered controller is:');
        disp(K);
    end
end