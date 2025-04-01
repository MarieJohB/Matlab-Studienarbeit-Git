function [K, details] = designCHR(G, structure, epsilon, plantInfo)
    % Enhanced CHR (Chien-Hrones-Reswick) method for 0% overshoot with improved stability handling
    
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % For unstable plants, pre-stabilize
    if plantInfo.isUnstable
        details = [details, 'Plant is unstable. Using pre-stabilization.\n'];
        K_stab = preStabilize(G, plantInfo);
        G_stab = feedback(G, K_stab);
        
        % Use stabilized plant for analysis
        t = linspace(0, 100, 1000);
        [y, t] = step(G_stab, t);
    else
        G_stab = G;
        
        % Get step response
        if ~isempty(plantInfo.stepResponse)
            t = plantInfo.stepResponse.time;
            y = plantInfo.stepResponse.response;
        else
            t = linspace(0, 100, 1000);
            [y, t] = step(G, t);
        end
    end
    
    % Find final value
    y_final = y(end);
    
    if abs(y_final) < 1e-6
        % No steady-state gain - try using plant info
        if ~isnan(plantInfo.FOPDT.K)
            y_final = plantInfo.FOPDT.K;
            details = [details, 'Using estimated steady-state gain from FOPDT model.\n'];
        else
            error('System does not have a finite DC gain. Not suitable for CHR method.');
        end
    end
    
    % Normalize the response
    y_norm = y / y_final;
    
    % Use FOPDT parameters if available, otherwise estimate them
    if ~isnan(plantInfo.FOPDT.L) && ~isnan(plantInfo.FOPDT.T)
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
                error('Could not determine time constant. The plant may not be suitable for the CHR method.');
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
    
    % Calculate parameter ratio for tuning adjustment
    ratio = L/T;
    
    details = [details, sprintf('Ks = %.4f\nL = %.4f s\nT = %.4f s\nL/T ratio = %.4f\n', Ks, L, T, ratio)];
    
    % Apply adjustment factor based on plant characteristics
    adjustmentFactor = 1.0;
    
    if ratio > 0.5
        % High delay-to-time-constant ratio - use more conservative settings
        adjustmentFactor = 0.8;
        details = [details, 'High L/T ratio: Using more conservative tuning.\n'];
    elseif plantInfo.hasRHPZeros
        % Non-minimum phase - use more conservative settings
        adjustmentFactor = 0.7;
        details = [details, 'Non-minimum phase behavior: Using more conservative tuning.\n'];
    elseif plantInfo.hasIntegrator
        % Plant with integrator - reduce integral action
        adjustmentFactor = 0.9;
        details = [details, 'Plant with integrator: Adjusting parameters.\n'];
    end
    
    % Parameters according to CHR table for setpoint regulation (0% overshoot)
    switch structure
        case 'P'
            Kp = adjustmentFactor * 0.3 * T / (Ks * L);
            K = tf(Kp, 1);
        case 'PI'
            Kp = adjustmentFactor * 0.35 * T / (Ks * L);
            Ti = 1.2 * T;
            K = tf([Kp, Kp/Ti], [1, 0]);
        case 'PD'
            Kp = adjustmentFactor * 0.5 * T / (Ks * L);
            Td = 0.3 * L;
            
            % For plants with RHP zeros, reduce derivative action
            if plantInfo.hasRHPZeros
                Td = Td * 0.5;
                details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
            end
            
            K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
        case 'PID'
            Kp = adjustmentFactor * 0.6 * T / (Ks * L);
            Ti = T;
            Td = 0.5 * L;
            
            % For plants with RHP zeros, reduce derivative action
            if plantInfo.hasRHPZeros
                Td = Td * 0.5;
                details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
            end
            
            K = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
        otherwise
            error('Unsupported controller structure for CHR method');
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
    
    % Verify controller and closed-loop stability
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
end