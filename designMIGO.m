function [K, details] = designMIGO(G, structure, robustness, epsilon, plantInfo)
    % Enhanced MIGO (M-constrained Integral Gain Optimization) method with better handling of unstable systems
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % Extract needed parameters
    epsilon = options.epsilon;
    robustness = options.robustness;
    
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
    
    % For unstable plants, first pre-stabilize
    if plantInfo.isUnstable
        details = [details, 'Plant is unstable. Using pre-stabilization before MIGO method.\n'];
        
        % Create a robust stabilizing controller for the unstable plant
        K_stab = designRobustStabilizingController(G, plantInfo);
        
        % Create a stabilized version of the plant
        G_stab = feedback(G * K_stab, 1);
        G_for_design = G_stab;
        isPreStabilized = true;
    else
        G_for_design = G;
        isPreStabilized = false;
    end
    
    % Implementation based on controller structure
    try
        switch structure
            case 'P'
                % For P-controller: find maximum Kp under M_s constraint
                Kp = findOptimalKp(G_for_design, M_s);
                K_migo = tf(Kp, 1);
                details = [details, sprintf('Optimal P-controller gain: Kp = %.4f\n', Kp)];
                
            case 'PI'
                % For PI-controller: optimize Kp and Ki under M_s constraint
                [Kp, Ki] = findOptimalPI(G_for_design, M_s, plantInfo);
                K_migo = tf([Kp, Ki], [1, 0]);
                Ti = Kp/Ki;
                details = [details, sprintf('Optimal PI parameters:\nKp = %.4f\nKi = %.4f\nTi = %.4f\n', Kp, Ki, Ti)];
                
            case 'PD'
                % For PD-controller: optimize Kp and Kd under M_s constraint
                [Kp, Kd] = findOptimalPD(G_for_design, M_s, epsilon, plantInfo);
                Td = Kd/Kp;
                K_migo = tf([Kd, Kp], [epsilon*Td, 1]);
                details = [details, sprintf('Optimal PD parameters:\nKp = %.4f\nKd = %.4f\nTd = %.4f\n', Kp, Kd, Td)];
                
            case 'PID'
                % For PID-controller: optimize Kp, Ki, and Kd under M_s constraint
                [Kp, Ki, Kd] = findOptimalPID(G_for_design, M_s, epsilon, plantInfo);
                Ti = Kp/Ki;
                Td = Kd/Kp;
                K_migo = tf([Kd, Kp, Ki], [epsilon*Td, 1, 0]);
                details = [details, sprintf('Optimal PID parameters:\nKp = %.4f\nKi = %.4f\nKd = %.4f\nTi = %.4f\nTd = %.4f\n', ...
                          Kp, Ki, Kd, Ti, Td)];
                
            otherwise
                error('MIGO is only implemented for P, PI, PD, and PID controllers');
        end
        
        % For pre-stabilized systems, combine with the stabilizing controller
        if isPreStabilized
            K_combined = series(K_migo, K_stab);
            
            try
                % Simplify the combined controller if possible
                K_combined = minreal(K_combined, 0.01);
                details = [details, 'Successfully simplified the combined controller.\n'];
            catch
                details = [details, 'Could not simplify the combined controller.\n'];
            end
            
            K = K_combined;
            details = [details, 'Combined with pre-stabilizing controller for final result.\n'];
        else
            K = K_migo;
        end
        
    catch ME
        details = [details, sprintf('Error in MIGO method: %s\n', ME.message)];
        
        % Fallback to conservative controller
        if isPreStabilized
            % If pre-stabilization worked, just use that controller
            K = K_stab;
            details = [details, 'Using pre-stabilizing controller as fallback.\n'];
        else
            % Create a conservative controller based on structure
            switch structure
                case 'P'
                    K = tf(0.1, 1);
                case 'PI'
                    K = tf([0.1, 0.01], [1, 0]);
                case 'PD'
                    K = tf([0.02, 0.1], [epsilon, 1]);
                case 'PID'
                    K = tf([0.02, 0.1, 0.01], [epsilon, 1, 0]);
                otherwise
                    K = tf(0.1, 1);
            end
            details = [details, 'Using default conservative controller as fallback.\n'];
        end
    end
    
    % Verify controller stability
    try
        K_poles = pole(K);
        if any(real(K_poles) > 0)
            details = [details, 'WARNING: Controller has unstable poles. Applying stabilization.\n'];
            
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
end