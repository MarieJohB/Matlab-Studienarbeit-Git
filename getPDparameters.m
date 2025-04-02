function [Kp, Kd] = getPDparameters(K)
    % Extract PD parameters from a general transfer function controller
    s = tf('s');
    [num, den] = tfdata(K, 'v');
    
    if length(den) == 1 && den(1) ~= 0
        % Proper transfer function
        if length(num) >= 2
            % PD controller form
            Kp = num(2) / den(1);
            Kd = num(1) / den(1);
        else
            % Only proportional term
            Kp = num(1) / den(1);
            Kd = 0;
        end
    else
        % Approximate by evaluating at specific frequencies
        w1 = 0.1;   % Lower frequency for proportional effect
        w2 = 10;    % Higher frequency for derivative effect
        
        G_w1 = abs(evalfr(K, 1j*w1));
        G_w2 = abs(evalfr(K, 1j*w2));
        
        Kp = G_w1;
        Kd = (G_w2 - G_w1) / (w2 - w1);
    end
end