function [K, details] = designZieglerNicholsStep(G, structure, epsilon, plantInfo)
    % Enhanced Ziegler-Nichols step method with adaptations for different plant types
    
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % For unstable plants, this method is not directly applicable
    if plantInfo.isUnstable
        % For unstable plants, use a different approach
        details = [details, 'Plant is unstable. Using modified approach.\n'];
        
        % Pre-stabilize the plant
        K_stab = preStabilize(G, plantInfo);
        G_stab = feedback(G, K_stab);
        
        % Now perform step response analysis on the stabilized plant
        t = linspace(0, 100, 1000);
        [y, t] = step(G_stab, t);
    else
        % For stable plants, use the standard approach
        G_stab = G;
        
        % Calculate step response
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
            error('System does not have a finite DC gain. Not suitable for Ziegler-Nichols step method.');
        end
    end
    
    % Normalize the response
    y_norm = y / y_final;
    
    % Process identification approach:
    % If we already have FOPDT parameters from plant analysis, use them
    if ~isnan(plantInfo.FOPDT.L) && ~isnan(plantInfo.FOPDT.T)
        L = plantInfo.FOPDT.L;
        T = plantInfo.FOPDT.T;
        Ks = plantInfo.FOPDT.K;
        
        details = [details, 'Using FOPDT parameters from plant analysis.\n'];
    else
        % Standard tangent method
        try
            % Find the derivative of the normalized response
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
            % Fallback to 63% method if tangent method fails
            try
                idx_63 = find(y_norm >= 0.632, 1);
                if isempty(idx_63)
                    idx_63 = find(y_norm >= 0.5, 1);
                    details = [details, 'Using 50% point instead of 63.2% due to limited response data.\n'];
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
            catch
                % If all methods fail, use conservative estimates
                L = 0.1 * t(end);
                T = 0.5 * t(end);
                Ks = y_final;
                
                details = [details, 'Using conservative estimates due to identification failures.\n'];
            end
        end
    end
    
    % Make sure parameters are reasonable
    L = max(L, 0.01);
    T = max(T, 0.1);
    
    % Calculate parameter ratio for tuning adjustment
    ratio = L/T;
    
    details = [details, sprintf('Ks = %.4f\nL = %.4f s\nT = %.4f s\nL/T ratio = %.4f\n', Ks, L, T, ratio)];
    
    % Modified parameters based on L/T ratio
    % For systems with significant delay, use more conservative settings
    if ratio > 0.5
        details = [details, 'High L/T ratio: Using more conservative tuning.\n'];
        conservativeFactor = 0.8;
    elseif ratio < 0.1
        details = [details, 'Low L/T ratio: Using more aggressive tuning.\n'];
        conservativeFactor = 1.2;
    else
        conservativeFactor = 1.0;
    end
    
    % Adjust for non-minimum phase behavior
    if plantInfo.hasRHPZeros
        conservativeFactor = conservativeFactor * 0.7;
        details = [details, 'Non-minimum phase behavior: Using more conservative tuning.\n'];
    end
    
    % Calculate parameters using adjusted Ziegler-Nichols rules
    switch structure
        case 'P'
            Kp = conservativeFactor * T / (Ks * L);
            K = tf(Kp, 1);
        case 'PI'
            Kp = conservativeFactor * 0.9 * T / (Ks * L);
            Ti = min(L / 0.3, 3 * L); % Limit Ti to avoid excessive integral action
            K = tf([Kp, Kp/Ti], [1, 0]);
        case 'PD'
            Kp = conservativeFactor * 1.2 * T / (Ks * L);
            Td = 0.5 * L;
            
            % For plants with RHP zeros, reduce derivative action
            if plantInfo.hasRHPZeros
                Td = Td * 0.5;
                details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
            end
            
            K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
        case 'PID'
            Kp = conservativeFactor * 1.2 * T / (Ks * L);
            Ti = 2 * L;
            Td = 0.5 * L;
            
            % For plants with RHP zeros, reduce derivative action
            if plantInfo.hasRHPZeros
                Td = Td * 0.5;
                details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
            end
            
            K = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
        otherwise
            error('Unsupported controller structure for Ziegler-Nichols step method');
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