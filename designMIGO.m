function [K, details] = designMIGO(G, structure, robustness, epsilon, plantInfo)
    % Enhanced MIGO (M-constrained Integral Gain Optimization) method
    
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % Set M-constraint based on robustness level
    switch robustness
        case 'Low'
            M_s = 2.0;  % Lower robustness, higher performance
        case 'Medium'
            M_s = 1.5;  % Medium robustness
        case 'High'
            M_s = 1.2;  % Higher robustness, lower performance
        otherwise
            M_s = 1.5;  % Default value
    end
    
    details = [details, sprintf('Robustness level: %s\nM_s constraint: %.2f\n', robustness, M_s)];
    
    % For unstable plants, pre-stabilize first
    if plantInfo.isUnstable
        details = [details, 'Plant is unstable. Using pre-stabilization approach.\n'];
        
        K_stab = preStabilize(G, plantInfo);
        G_stab = feedback(G, K_stab);
        
        % Continue MIGO design with stabilized plant
        G_for_design = G_stab;
    else
        G_for_design = G;
    end
    
    % Implementation based on controller structure
    switch structure
        case 'P'
            % For P-controller: find maximum Kp under M_s constraint
            Kp = findOptimalKp(G_for_design, M_s);
            K = tf(Kp, 1);
            details = [details, sprintf('Optimal P-controller gain: Kp = %.4f\n', Kp)];
            
        case 'PI'
            % For PI-controller: optimize Kp and Ki under M_s constraint
            [Kp, Ki] = findOptimalPI(G_for_design, M_s, plantInfo);
            K = tf([Kp, Ki], [1, 0]);
            Ti = Kp/Ki;
            details = [details, sprintf('Optimal PI parameters:\nKp = %.4f\nKi = %.4f\nTi = %.4f\n', Kp, Ki, Ti)];
            
        case 'PD'
            % For PD-controller: optimize Kp and Kd
            [Kp, Kd] = findOptimalPD(G_for_design, M_s, epsilon, plantInfo);
            K = tf([Kd, Kp], [epsilon*Kd, 1]);
            Td = Kd/Kp;
            details = [details, sprintf('Optimal PD parameters:\nKp = %.4f\nKd = %.4f\nTd = %.4f\n', Kp, Kd, Td)];
            
        case 'PID'
            % For PID-controller: optimize Kp, Ki, and Kd under M_s constraint
            [Kp, Ki, Kd] = findOptimalPID(G_for_design, M_s, epsilon, plantInfo);
            K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
            Ti = Kp/Ki;
            Td = Kd/Kp;
            details = [details, sprintf('Optimal PID parameters:\nKp = %.4f\nKi = %.4f\nKd = %.4f\nTi = %.4f\nTd = %.4f\n', ...
                      Kp, Ki, Kd, Ti, Td)];
            
        otherwise
            error('MIGO is only implemented for P, PI, PD, and PID controllers');
    end
    
    % For unstable plants, combine with pre-stabilizing controller
    if plantInfo.isUnstable
        K_combined = series(K, K_stab);
        
        % Simplify the combined controller if possible
        try
            K_combined = minreal(K_combined);
        catch
            % If simplification fails, use the original combination
        end
        
        K = K_combined;
        details = [details, 'Combined with pre-stabilizing controller for unstable plant.\n'];
    end
    
    % Verify controller stability
    try
        K_poles = pole(K);
        if any(real(K_poles) > 0)
            details = [details, 'Warning: Controller has unstable poles. Applying stabilization.\n'];
            
            % Force controller stability
            [num, den] = tfdata(K, 'v');
            p = roots(den);
            for i = 1:length(p)
                if real(p(i)) > 0
                    p(i) = -abs(real(p(i))) + imag(p(i))*1i;
                end
            end
            den_stable = poly(p);
            K = tf(num, den_stable);
        end
    catch
        % If pole analysis fails, continue with the designed controller
    end
    
    % Verify closed-loop stability
    try
        T = feedback(G * K, 1);
        cl_poles = pole(T);
        
        if any(real(cl_poles) > 0)
            details = [details, 'Warning: Closed-loop system is unstable. Adjusting controller.\n'];
            
            % Try to stabilize by reducing gain
            [num, den] = tfdata(K, 'v');
            K_adjusted = tf(num * 0.5, den);
            
            % Check if modification helps
            T_adj = feedback(G * K_adjusted, 1);
            
            if all(real(pole(T_adj)) < 0)
                K = K_adjusted;
                details = [details, 'Controller gain reduced by 50% to achieve stability.\n'];
            else
                % Try more aggressive reduction
                K_adjusted = tf(num * 0.2, den);
                T_adj = feedback(G * K_adjusted, 1);
                
                if all(real(pole(T_adj)) < 0)
                    K = K_adjusted;
                    details = [details, 'Controller gain reduced by 80% to achieve stability.\n'];
                else
                    details = [details, 'Could not stabilize system by gain reduction. Consider a different method.\n'];
                end
            end
        end
    catch
        details = [details, 'Could not verify closed-loop stability.\n'];
    end
end