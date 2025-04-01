function [K, details] = designAstrom(G, structure, epsilon, plantInfo)
    % Enhanced Åström method with improved parameter selection and stability handling
    
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
            error('System does not have a finite DC gain. Not suitable for Åström method.');
        end
    end
    
    % Static gain
    Ks = y_final;
    
    % Use FOPDT parameters if available, otherwise estimate them
    if ~isnan(plantInfo.FOPDT.L) && ~isnan(plantInfo.FOPDT.T)
        L = plantInfo.FOPDT.L;
        T = plantInfo.FOPDT.T;
    else
        % Estimate time constant and delay
        y_norm = y / y_final;
        
        % Find 63.2% point for T
        idx_63 = find(y_norm >= 0.632, 1);
        
        if isempty(idx_63)
            error('Could not determine time constant. The plant may not be suitable for the Åström method.');
        end
        
        % Estimate delay by comparing with first-order model plus delay
        idx_10 = find(y_norm >= 0.1, 1);
        if isempty(idx_10)
            L = 0.1;  % Default value
        else
            L = t(idx_10);  % Estimate delay as time at 10% rise
        end
        
        T = t(idx_63) - L;  % Time constant (63.2% time minus delay)
    end
    
    % Calculate parameter ratio for tuning adjustment
    ratio = L/T;
    
    details = [details, sprintf('Ks = %.4f\nL = %.4f s\nT = %.4f s\nL/T ratio = %.4f\n', Ks, L, T, ratio)];
    
    % Calculate controller parameters with improved Åström rules
    switch structure
        case 'P'
            if L < 0.5*T
                Kp = 0.3 * T / (Ks * L);
            elseif L < 2*T
                Kp = 0.2 * T / (Ks * L); 
            else
                Kp = 0.15 / Ks;
            end
            K = tf(Kp, 1);
        case 'PI'
            if L < 0.5*T
                Kp = 0.3 * T / (Ks * L);
            elseif L < 2*T
                Kp = 0.25 * T / (Ks * L); 
            else
                Kp = 0.15 / Ks;
            end
            
            if L < 0.1*T
                Ti = 8 * L;
            elseif L < 2*T
                Ti = 0.8 * T;
            else
                Ti = 0.4 * (L + T);
            end
            
            K = tf([Kp, Kp/Ti], [1, 0]);
        case 'PD'
            Kp = 0.4 * T / (Ks * L);
            Td = 0.4 * L;
            
            % For plants with RHP zeros, reduce derivative action
            if plantInfo.hasRHPZeros
                Td = Td * 0.5;
                details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
            end
            
            K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
        case 'PID'
            if plantInfo.hasRHPZeros
                details = [details, 'Adjusting PID for non-minimum phase plant.\n'];
                Kp = 0.25 * T / (Ks * L);
                Ti = 1.2 * T;
                Td = 0.2 * L;
            elseif L < 0.1*T
                Kp = 0.3 * T / (Ks * L);
                Ti = 8 * L;
                Td = 0.25 * L;
            elseif L < 2*T
                Kp = 0.3 * T / (Ks * L);
                Ti = 0.8 * T;
                Td = 0.2 * L;
            else
                Kp = 0.15 / Ks;
                Ti = 0.4 * (L + T);
                Td = 0.15 * L;
            end
            
            % For plants with integrator, adjust integral term
            if plantInfo.hasIntegrator
                details = [details, 'Adjusting for integrator in plant.\n'];
                Ti = Ti * 2;
            end
            
            K = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
        otherwise
            error('Unsupported controller structure for Åström method');
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