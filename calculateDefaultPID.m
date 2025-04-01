function defaultPID = calculateDefaultPID(plantInfo)
    % CALCULATEDEFAULTPID Calculate default PID parameters based on plant characteristics
    
    defaultPID = struct();
    
    if ~isnan(plantInfo.FOPDT.K) && ~isnan(plantInfo.FOPDT.T) && ~isnan(plantInfo.FOPDT.L)
        % Use FOPDT parameters
        K = plantInfo.FOPDT.K;
        T = plantInfo.FOPDT.T;
        L = plantInfo.FOPDT.L;
        
        if plantInfo.isUnstable
            % For unstable plants
            defaultPID.Kp = plantInfo.stabilizingGain * 1.2;
            defaultPID.Ki = defaultPID.Kp * 0.1;
            defaultPID.Kd = defaultPID.Kp * 0.1;
        else
            % Standard SIMC tuning rules (Skogestad)
            c = 4; % Tuning factor (4 for normal, higher for more robustness)
            
            % PI parameters
            defaultPID.Kp = (1/K) * (T / (c*L));
            defaultPID.Ti = min(T, 4*L);
            defaultPID.Ki = defaultPID.Kp / defaultPID.Ti;
            
            % D parameter
            if L > 0.1*T
                defaultPID.Td = L;
            else
                defaultPID.Td = 0;
            end
            defaultPID.Kd = defaultPID.Kp * defaultPID.Td;
            
            % Filter for derivative term
            defaultPID.N = 10;
        end
    else
        % Fallback to rule-of-thumb parameters
        try
            if ~plantInfo.isUnstable
                % Get gain at a frequency near the bandwidth
                if ~isnan(plantInfo.bandwidth)
                    w = plantInfo.bandwidth;
                else
                    % Try to find a reasonable frequency
                    p = plantInfo.poles;
                    stable_poles = p(real(p) < 0);
                    if ~isempty(stable_poles)
                        w = min(abs(real(stable_poles)));
                    else
                        w = 1.0; % Default
                    end
                end
                
                % Calculate gain at this frequency
                try
                    g = abs(evalfr(G, 1j*w));
                    defaultPID.Kp = 0.5 / g;
                catch
                    defaultPID.Kp = 1.0;
                end
                
                defaultPID.Ti = 1 / w;
                defaultPID.Ki = defaultPID.Kp / defaultPID.Ti;
                
                defaultPID.Td = 0.25 / w;
                defaultPID.Kd = defaultPID.Kp * defaultPID.Td;
                
                defaultPID.N = 10;
            else
                % For unstable plants
                defaultPID.Kp = plantInfo.stabilizingGain * 1.2;
                defaultPID.Ki = defaultPID.Kp * 0.1;
                defaultPID.Kd = defaultPID.Kp * 0.1;
                defaultPID.Ti = defaultPID.Kp / defaultPID.Ki;
                defaultPID.Td = defaultPID.Kd / defaultPID.Kp;
                defaultPID.N = 10;
            end
        catch
            % Default values if all else fails
            defaultPID.Kp = 1.0;
            defaultPID.Ki = 0.1;
            defaultPID.Kd = 0.1;
            defaultPID.Ti = 10;
            defaultPID.Td = 0.1;
            defaultPID.N = 10;
        end
    end
end