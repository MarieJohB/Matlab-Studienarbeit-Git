function Kp = findOptimalKp(G, M_s)
    % Find the maximum P-controller gain that satisfies the M_s constraint
    
    % Define search range
    Kp_min = 0.01;
    Kp_max = 100;
    Kp_values = logspace(log10(Kp_min), log10(Kp_max), 50);
    
    Kp_best = Kp_min;
    
    for i = 1:length(Kp_values)
        Kp = Kp_values(i);
        K_test = tf(Kp, 1);
        L = G * K_test;
        S = feedback(1, L);
        
        try
            % Calculate maximum sensitivity
            [peakgain, ~, ~] = getPeakGain(S);
            
            if peakgain <= M_s && Kp > Kp_best
                Kp_best = Kp;
            end
        catch
            % Skip if peak gain calculation fails
            continue;
        end
    end
    
    Kp = Kp_best;
end