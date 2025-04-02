function [Kp, Ki] = getPIparameters(K)
    % Extract PI parameters from a general transfer function controller
    s = tf('s');
    [num, den] = tfdata(K, 'v');
    
    if length(den) == 1 && den(1) ~= 0
        % Proper transfer function
        if length(num) >= 2
            if den(1) == 0
                % Pure integral term exists
                Kp = num(1);
                Ki = num(2);
            else
                % Standard case
                Kp = num(1) / den(1);
                Ki = (num(length(num)) / den(1)) * den(length(den));
            end
        else
            % Only proportional term
            Kp = num(1) / den(1);
            Ki = Kp / 10;  % Default integral action
        end
    else
        % Approximate by evaluating at specific frequencies
        w1 = 0.01;  % Low frequency for integral effect
        w2 = 10;    % Higher frequency for proportional effect
        
        G_w1 = abs(evalfr(K, 1j*w1));
        G_w2 = abs(evalfr(K, 1j*w2));
        
        Kp = G_w2;
        Ki = G_w1 * w1;
    end
end