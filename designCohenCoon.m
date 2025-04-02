function [K, details] = designCohenCoon(G, structure, options, plantInfo)
    % Enhanced Cohen-Coon method with pre-stabilization for unstable plants
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % Extract needed parameters
    epsilon = options.epsilon;
    
    % For unstable systems, first pre-stabilize
    if plantInfo.isUnstable
        details = [details, 'Plant is unstable. Using pre-stabilization before Cohen-Coon method.\n'];
        
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
    
    % Apply the Cohen-Coon method
    try
        % Get or generate step response
        if ~isempty(plantInfo.stepResponse) && ~isPreStabilized
            t = plantInfo.stepResponse.time;
            y = plantInfo.stepResponse.response;
        else
            t = linspace(0, 100, 1000);
            [y, t] = step(G_for_design, t);
            
            % Reshape y if it's multidimensional
            if size(y, 2) > 1
                y = y(:, 1, 1);
            end
        end
        
        % Make sure we have a valid step response
        y_final = y(end);
        if abs(y_final) < 1e-6
            error('System does not have a finite DC gain');
        end
        
        % Normalize the response
        y_norm = y / y_final;
        
        % Determine FOPDT parameters
        if ~isnan(plantInfo.FOPDT.L) && ~isnan(plantInfo.FOPDT.T) && ~isPreStabilized
            % Use FOPDT parameters if already available
            L = plantInfo.FOPDT.L;
            T = plantInfo.FOPDT.T;
            Ks = plantInfo.FOPDT.K;
        else
            try
                % Determine the tangent at the inflection point
                dy = diff(y_norm) ./ diff(t);
                
                % Find the point with maximum slope (inflection point)
                [max_slope, idx_max_slope] = max(dy);
                
                % Calculate parameters of the tangent: y = m*t + b
                m = max_slope;
                b = y_norm(idx_max_slope) - m * t(idx_max_slope);
                
                % Intersections with y=0 and y=1
                t_0 = -b / m;                 % Intersection with y=0 (t-axis)
                t_1 = (1 - b) / m;            % Intersection with y=1
                
                % Determine dead time L and rise time T
                L = t_0;
                T = t_1 - t_0;
                
                % Static gain
                Ks = y_final;
            catch
                % Fallback method if tangent method fails
                idx_63 = find(y_norm >= 0.632, 1);
                if isempty(idx_63)
                    error('Could not determine time constant. The plant may not be suitable for the Cohen-Coon method.');
                end
                
                idx_10 = find(y_norm >= 0.1, 1);
                if isempty(idx_10)
                    L = 0.1 * t(idx_63);
                else
                    L = t(idx_10);
                end
                
                T = t(idx_63) - L;
                Ks = y_final;
                
                details = [details, 'Using 63.2% method due to tangent method failure.\n'];
            end
        end
        
        % Parameter tau = T/L, helpful for Cohen-Coon formulas
        tau = T/L;
        
        details = [details, sprintf('Ks = %.4f\nL = %.4f s\nT = %.4f s\ntau = %.4f\n', Ks, L, T, tau)];
        
        % Apply adjustment factor based on plant characteristics
        adjustmentFactor = 1.0;
        
        if tau < 3
            % Low tau - reduce controller gain for stability
            adjustmentFactor = 0.8;
            details = [details, 'Low tau value: Reducing controller gain for stability.\n'];
        elseif plantInfo.hasRHPZeros
            % Non-minimum phase - use more conservative settings
            adjustmentFactor = 0.7;
            details = [details, 'Non-minimum phase behavior: Using more conservative tuning.\n'];
        elseif plantInfo.hasDelay
            adjustmentFactor = 0.85;
            details = [details, 'Detected significant delay: Adjusting parameters.\n'];
        elseif plantInfo.isHighOrder
            adjustmentFactor = 0.9;
            details = [details, 'High-order system: Adjusting for robustness.\n'];
        end
        
        % Improved Cohen-Coon formulas with adjustments
        switch structure
            case 'P'
                Kp = adjustmentFactor * (1/Ks) * (1 + (1/3)*tau);
                K_cc = tf(Kp, 1);
                details = [details, sprintf('P controller: Kp = %.4f\n', Kp)];
                
            case 'PI'
                Kp = adjustmentFactor * (1/Ks) * (0.9 + (1/12)*tau);
                Ti = L * ((30 + 3*tau)/(9 + 20*tau));
                Ki = Kp / Ti;
                
                % For plants with integrator, adjust Ti
                if plantInfo.hasIntegrator
                    Ti = Ti * 1.5;
                    Ki = Kp / Ti;
                    details = [details, 'Adjusted Ti for plant with integrator.\n'];
                end
                
                K_cc = tf([Kp, Ki], [1, 0]);
                details = [details, sprintf('PI controller: Kp = %.4f, Ti = %.4f, Ki = %.4f\n', Kp, Ti, Ki)];
                
            case 'PD'
                Kp = adjustmentFactor * (1/Ks) * (1.25 * (1 + (1/6)*tau));
                Td = L * ((6 - 2*tau)/(22 + 3*tau));
                
                % For plants with RHP zeros, reduce derivative action
                if plantInfo.hasRHPZeros
                    Td = Td * 0.5;
                    details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
                end
                
                % Ensure Td is positive
                Td = max(Td, 0.05*L);
                Kd = Kp * Td;
                
                K_cc = tf([Kd, Kp], [epsilon*Td, 1]);
                details = [details, sprintf('PD controller: Kp = %.4f, Td = %.4f, Kd = %.4f\n', Kp, Td, Kd)];
                
            case 'PID'
                Kp = adjustmentFactor * (1/Ks) * (1.35 + (1/4)*tau);
                Ti = L * ((32 + 6*tau)/(13 + 8*tau));
                Td = L * (4/(11 + 2*tau));
                
                % For plants with RHP zeros, reduce derivative action
                if plantInfo.hasRHPZeros
                    Td = Td * 0.5;
                    details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
                end
                
                % For plants with integrator, adjust Ti
                if plantInfo.hasIntegrator
                    Ti = Ti * 1.5;
                    details = [details, 'Adjusted Ti for plant with integrator.\n'];
                end
                
                Ki = Kp / Ti;
                Kd = Kp * Td;
                
                K_cc = tf([Kd, Kp, Ki], [epsilon*Td, 1, 0]);
                details = [details, sprintf('PID controller: Kp = %.4f, Ti = %.4f, Td = %.4f, Ki = %.4f, Kd = %.4f\n', ...
                          Kp, Ti, Td, Ki, Kd)];
                
            otherwise
                error('Unsupported controller structure for Cohen-Coon method');
        end
        
        % For pre-stabilized systems, combine with the stabilizing controller
        if isPreStabilized
            K_combined = series(K_cc, K_stab);
            
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
            K = K_cc;
        end
        
    catch ME
        details = [details, sprintf('Error in Cohen-Coon method: %s\n', ME.message)];
        
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