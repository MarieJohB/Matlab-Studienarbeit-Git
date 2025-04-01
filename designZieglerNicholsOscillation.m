function [K, details] = designZieglerNicholsOscillation(G, structure, epsilon, plantInfo)
    % Improved Ziegler-Nichols oscillation method with pre-stabilization for unstable plants
    
    % Add plant information to details
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % For unstable plants, pre-stabilize before Z-N tuning
    if plantInfo.isUnstable
        K_stab = preStabilize(G, plantInfo);
        details = [details, sprintf('Plant is unstable. Using pre-stabilization.\n')];
        G_stab = feedback(G, K_stab);
    else
        G_stab = G;
    end
    
    % Determine critical gain and period
    k_krit = 0;
    T_krit = 0;
    
    try
        % Use frequency domain approach for critical point
        w = logspace(-2, 4, 1000);
        [mag, phase] = bode(G_stab, w);
        mag = squeeze(mag);
        phase = squeeze(phase);
        
        % Find where phase crosses -180 degrees
        phase_crossings = [];
        for i = 1:length(phase)-1
            if (phase(i) > -180 && phase(i+1) < -180) || (phase(i) < -180 && phase(i+1) > -180)
                phase_crossings(end+1) = i;
            end
        end
        
        if ~isempty(phase_crossings)
            % Linear interpolation to find precise crossing frequency
            for crossing = phase_crossings
                w1 = w(crossing);
                w2 = w(crossing+1);
                phase1 = phase(crossing);
                phase2 = phase(crossing+1);
                
                w_cross = w1 + (w2 - w1) * (-180 - phase1) / (phase2 - phase1);
                
                % Get magnitude at this frequency using interpolation
                mag1 = mag(crossing);
                mag2 = mag(crossing+1);
                mag_cross = mag1 + (mag2 - mag1) * (w_cross - w1) / (w2 - w1);
                
                % Critical gain is reciprocal of magnitude at -180 degrees
                k_krit_candidate = 1 / mag_cross;
                T_krit_candidate = 2*pi / w_cross;
                
                % Use the first valid crossing point
                if k_krit == 0 || k_krit_candidate < k_krit
                    k_krit = k_krit_candidate;
                    T_krit = T_krit_candidate;
                end
            end
        end
        
        % If frequency domain approach failed, use time-domain approach
        if k_krit == 0
            % Start with a small gain and increase it iteratively
            found = false;
            step_size = 0.1;
            max_iterations = 2000;
            
            for iteration = 1:max_iterations
                k = iteration * step_size;
                
                % Test stability with current k
                closed_loop = feedback(G_stab*k, 1);
                poles = pole(closed_loop);
                
                % Check for poles on the imaginary axis
                realParts = real(poles);
                imagParts = imag(poles);
                
                close_to_imag_axis = find(abs(realParts) < 0.001 & imagParts ~= 0);
                
                if ~isempty(close_to_imag_axis)
                    % Found poles near stability boundary
                    k_krit = k;
                    
                    % Calculate period
                    idx = find(abs(realParts) < 0.001 & imagParts > 0, 1);
                    if ~isempty(idx)
                        omega = imagParts(idx);
                        T_krit = 2*pi/omega;
                        found = true;
                        break;
                    end
                elseif any(realParts > 0)
                    % System became unstable
                    
                    % Reduce k incrementally to find boundary
                    for back_step = 1:10
                        k_test = k - back_step * (step_size/10);
                        
                        closed_loop = feedback(G_stab*k_test, 1);
                        poles = pole(closed_loop);
                        realParts = real(poles);
                        imagParts = imag(poles);
                        
                        if all(realParts < 0) && any(abs(realParts) < 0.01 & imagParts ~= 0)
                            k_krit = k_test;
                            
                            idx = find(abs(realParts) < 0.01 & imagParts > 0, 1);
                            if ~isempty(idx)
                                omega = imagParts(idx);
                                T_krit = 2*pi/omega;
                                found = true;
                                break;
                            end
                        end
                    end
                    
                    if found
                        break;
                    else
                        % Assume we've passed the critical point
                        k_krit = k - step_size;
                        
                        closed_loop = feedback(G_stab*k_krit, 1);
                        poles = pole(closed_loop);
                        imagParts = imag(poles(abs(real(poles)) < 0.1));
                        
                        if ~isempty(imagParts) && any(imagParts > 0)
                            omega = max(imagParts(imagParts > 0));
                            T_krit = 2*pi/omega;
                            found = true;
                            break;
                        end
                    end
                    
                    break;
                end
            end
        end
        
        if k_krit == 0
            % Fallback estimation based on plant characteristics
            if ~isnan(plantInfo.FOPDT.L) && ~isnan(plantInfo.FOPDT.T)
                % Approximate critical values based on FOPDT model
                k_krit = (plantInfo.FOPDT.T / plantInfo.FOPDT.K / plantInfo.FOPDT.L) * pi / 2;
                T_krit = 4 * plantInfo.FOPDT.L;
                details = [details, sprintf('Using estimated critical parameters from FOPDT model.\n')];
            else
                error('Could not determine critical gain. The plant may not be suitable for the Ziegler-Nichols oscillation method.');
            end
        end
        
    catch ME
        % If all else fails, use approximation based on poles
        if plantInfo.isHighOrder
            % For higher-order systems, use conservative approximation
            k_krit = 2.0;
            T_krit = 1.0;
            details = [details, sprintf('Critical parameters estimated due to error: %s\n', ME.message)];
        else
            error('Error applying Ziegler-Nichols oscillation method: %s', ME.message);
        end
    end
    
    % Controller parameters according to Ziegler-Nichols table
    details = [details, sprintf('k_krit = %.4f\nT_krit = %.4f s\n', k_krit, T_krit)];
    
    switch structure
        case 'P'
            Kp = 0.5 * k_krit;
            K = tf(Kp, 1);
        case 'PI'
            Kp = 0.45 * k_krit;
            Ti = 0.85 * T_krit;
            K = tf([Kp, Kp/Ti], [1, 0]);
        case 'PD'
            Kp = 0.5 * k_krit;
            Td = 0.12 * T_krit;
            K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
        case 'PID'
            % Use modified parameters for better robustness
            Kp = 0.6 * k_krit;
            Ti = 0.5 * T_krit;
            Td = 0.125 * T_krit;
            
            % For plants with RHP zeros, reduce derivative action
            if plantInfo.hasRHPZeros
                Td = Td * 0.5;
                details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
            end
            
            K = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
        otherwise
            error('Unsupported controller structure for Ziegler-Nichols oscillation method');
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
            
            % Force controller stability by modifying unstable poles
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