% Simple Ziegler-Nichols tuning function
function [Kp, Ki, Kd] = simpleZieglerNichols(sys)
    % This is a simplified version of Z-N tuning
    % In a full implementation, you'd identify Kcr and Pcr
    
    % Try to estimate ultimate gain and period
    try
        % Get frequency response
        w = logspace(-2, 2, 1000);
        [mag, phase] = bode(sys, w);
        
        % Find frequency where phase is -180 degrees
        phase_180_idx = find(phase(:) <= -180, 1);
        
        if ~isempty(phase_180_idx)
            w_180 = w(phase_180_idx);
            mag_180 = mag(phase_180_idx);
            
            % Ultimate gain and period
            Ku = 1/mag_180;
            Pu = 2*pi/w_180;
            
            % PID Ziegler-Nichols tuning
            Kp = 0.6 * Ku;
            Ti = 0.5 * Pu;
            Td = 0.125 * Pu;
            
            Ki = Kp / Ti;
            Kd = Kp * Td;
        else
            % Default conservative values if no phase crossover
            Kp = 1;
            Ki = 0.1;
            Kd = 0.1;
        end
    catch
        % Fallback values
        Kp = 1;
        Ki = 0.1;
        Kd = 0.1;
    end
end