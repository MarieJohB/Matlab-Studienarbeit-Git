function [Kp, Ki, Kd] = getPIDparameters(K)
    % Extract PID parameters from a general transfer function controller
    s = tf('s');
    [num, den] = tfdata(K, 'v');
    
    if length(den) == 1 && den(1) ~= 0
        % Proper transfer function
        if length(num) >= 3
            % Full PID controller form
            Kd = num(1) / den(1);
            Kp = num(2) / den(1);
            Ki = num(3) / den(1);
        elseif length(num) == 2
            % PD or PI form
            Kp = num(1) / den(1);
            
            % Determine if it's PI or PD
            if num(1) == 0
                Kd = 0;
                Ki = num(2) / den(1);
            else
                Kd = num(1) / den(1);
                Ki = 0;
            end
        else
            % Only proportional term
            Kp = num(1) / den(1);
            Ki = 0;
            Kd = 0;
        end
    else
        % Approximate by evaluating at specific frequencies
        w1 = 0.01;  % Low frequency for integral effect
        w2 = 1;     % Mid frequency for proportional effect
        w3 = 100;   % High frequency for derivative effect
        
        G_w1 = abs(evalfr(K, 1j*w1));
        G_w2 = abs(evalfr(K, 1j*w2));
        G_w3 = abs(evalfr(K, 1j*w3));
        
        Ki = G_w1 * w1;
        Kp = G_w2;
        Kd = (G_w3 - G_w2) / (w3 - w2);
    end
end