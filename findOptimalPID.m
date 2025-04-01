function [Kp, Ki, Kd] = findOptimalPID(G, M_s, epsilon, plantInfo)
    % Find optimal PID parameters under M_s constraint
    % Multi-stage approach for PID design
    
    % 1. Find good PI parameters
    [Kp_pi, Ki_pi] = findOptimalPI(G, M_s * 1.1, plantInfo);
    
    % 2. Add derivative action starting from PI solution
    Td_values = logspace(-2, 0.5, 10);
    
    best_Kp = Kp_pi;
    best_Ki = Ki_pi;
    best_Td = 0;
    best_score = -Inf;
    
    for j = 1:length(Td_values)
        Td = Td_values(j);
        Kd = Kp_pi * Td;
        
        K_test = tf([Kd, Kp_pi, Ki_pi], [epsilon*Td, 1, 0]);
        L = G * K_test;
        
        try
            % Calculate closed-loop response
            S = feedback(1, L);
            [peakgain, ~, ~] = getPeakGain(S);
            [~, Pm] = margin(L);
            
            if peakgain <= M_s
                % Score based on phase margin and integral gain
                score = Ki_pi * sqrt(Kd) * (Pm/100);
                
                if score > best_score
                    best_score = score;
                    best_Td = Td;
                end
            end
        catch
            % Skip if analysis fails
            continue;
        end
    end
    
    % 3. Fine-tune the PID controller
    best_Kp = Kp_pi;
    best_Ki = Ki_pi;
    best_Kd = best_Kp * best_Td;
    
    % Adjust parameters based on plant characteristics
    if plantInfo.hasIntegrator
        % For plants with integrator, reduce integral action
        best_Ki = best_Ki * 0.7;
    end
    
    if plantInfo.hasRHPZeros
        % For non-minimum phase plants, reduce derivative action
        best_Kd = best_Kd * 0.7;
    end
    
    if plantInfo.isHighOrder
        % For high-order plants, add more filtering
        epsilon = epsilon * 1.5;
    end
    
    Kp = best_Kp;
    Ki = best_Ki;
    Kd = best_Kd;
end