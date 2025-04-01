function [Kp, Ki] = findOptimalPI(G, M_s, plantInfo)
    % Find optimal PI controller parameters under M_s constraint
    
    % Grid search approach with refinement
    % Start with coarse grid
    Kp_values = logspace(-2, 2, 15);
    Ti_values = logspace(-2, 2, 15);
    
    best_Kp = 0.1;
    best_Ti = 1.0;
    best_score = -Inf;
    
    % First pass: coarse grid search
    for i = 1:length(Kp_values)
        for j = 1:length(Ti_values)
            Kp = Kp_values(i);
            Ti = Ti_values(j);
            Ki = Kp / Ti;
            
            K_test = tf([Kp, Ki], [1, 0]);
            L = G * K_test;
            S = feedback(1, L);
            
            try
                % Calculate maximum sensitivity
                [peakgain, ~, ~] = getPeakGain(S);
                
                if peakgain <= M_s
                    % Score based on combined performance metrics
                    [~, Pm] = margin(L);
                    
                    % Calculate step response metrics if possible
                    try
                        T = feedback(G * K_test, 1);
                        step_info = stepinfo(T);
                        settling_time = step_info.SettlingTime;
                        overshoot = step_info.Overshoot;
                        
                        % Calculate score based on multiple criteria
                        % - Higher Ki for better tracking
                        % - Higher phase margin for stability
                        % - Lower overshoot and settling time for better performance
                        score = Ki * (Pm/100) * (100/(1 + overshoot)) * (10/(1 + settling_time));
                    catch
                        % If step response fails, use simpler scoring
                        score = Ki * (Pm/100);
                    end
                    
                    if score > best_score
                        best_score = score;
                        best_Kp = Kp;
                        best_Ti = Ti;
                    end
                end
            catch
                % Skip if analysis fails
                continue;
            end
        end
    end
    
    % Second pass: refined search around best point
    Kp_min = best_Kp * 0.5;
    Kp_max = best_Kp * 2.0;
    Ti_min = best_Ti * 0.5;
    Ti_max = best_Ti * 2.0;
    
    Kp_refined = linspace(Kp_min, Kp_max, 10);
    Ti_refined = linspace(Ti_min, Ti_max, 10);
    
    for i = 1:length(Kp_refined)
        for j = 1:length(Ti_refined)
            Kp = Kp_refined(i);
            Ti = Ti_refined(j);
            Ki = Kp / Ti;
            
            K_test = tf([Kp, Ki], [1, 0]);
            L = G * K_test;
            S = feedback(1, L);
            
            try
                % Calculate maximum sensitivity
                [peakgain, ~, ~] = getPeakGain(S);
                
                if peakgain <= M_s
                    % Score based on phase margin and integral gain
                    [~, Pm] = margin(L);
                    score = Ki * (Pm/100);
                    
                    if score > best_score
                        best_score = score;
                        best_Kp = Kp;
                        best_Ti = Ti;
                    end
                end
            catch
                % Skip if analysis fails
                continue;
            end
        end
    end
    
    % Adapt to plant characteristics
    if plantInfo.hasIntegrator
        % For plants with integrator, reduce integral action
        best_Ti = best_Ti * 2;
    end
    
    if plantInfo.hasRHPZeros
        % For non-minimum phase plants, reduce controller gain
        best_Kp = best_Kp * 0.8;
    end
    
    Kp = best_Kp;
    Ki = Kp / best_Ti;
end