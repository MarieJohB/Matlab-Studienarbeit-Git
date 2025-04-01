function [Kp, Kd] = findOptimalPD(G, M_s, epsilon, plantInfo)
    % Find optimal PD controller parameters under M_s constraint
    
    % Grid search approach
    Kp_values = logspace(-2, 2, 15);
    Td_values = logspace(-2, 1, 15);
    
    best_Kp = 0.1;
    best_Td = 0.1;
    best_Pm = 0;
    
    for i = 1:length(Kp_values)
        for j = 1:length(Td_values)
            Kp = Kp_values(i);
            Td = Td_values(j);
            Kd = Kp * Td;
            
            % Fix: Use epsilon*Td instead of epsilon*Kd in the denominator
            K_test = tf([Kd, Kp], [epsilon*Td, 1]);
            L = G * K_test;
            
            try
                % Calculate maximum sensitivity and phase margin
                S = feedback(1, L);
                [peakgain, ~, ~] = getPeakGain(S);
                [~, Pm] = margin(L);
                
                if peakgain <= M_s && Pm > best_Pm
                    best_Pm = Pm;
                    best_Kp = Kp;
                    best_Td = Td;
                end
            catch
                % Skip if analysis fails
                continue;
            end
        end
    end
    
    % Adapt to plant characteristics
    if plantInfo.hasRHPZeros
        % For non-minimum phase plants, reduce derivative action
        best_Td = best_Td * 0.7;
    end
    
    if plantInfo.isHighOrder
        % For high-order plants, add more filtering
        epsilon = epsilon * 1.5;
    end
    
    Kp = best_Kp;
    Kd = Kp * best_Td;
end