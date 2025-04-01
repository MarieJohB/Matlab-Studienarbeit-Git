function [K, details] = designIMC(G, structure, settlingTime, epsilon, plantInfo)
    % Enhanced IMC (Internal Model Control) method with better handling of stability and non-minimum phase
    
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    details = [details, sprintf('Target Settling Time: %.2f s\n', settlingTime)];
    
    % IMC filter time constant based on settling time
    % Estimate: For a first-order system, settling time is about 4*time constant
    lambda = settlingTime / 4;
    
    details = [details, sprintf('Lambda = %.4f (based on settling time)\n', lambda)];
    
    % For unstable plants, use a stabilizing inner loop first
    if plantInfo.isUnstable
        details = [details, 'Plant is unstable. Using pre-stabilization approach.\n'];
        
        % Apply pre-stabilization
        K_stab = preStabilize(G, plantInfo);
        G_stab = feedback(G, K_stab);
        
        % Continue IMC design with stabilized plant
        G_for_design = G_stab;
    else
        G_for_design = G;
    end
    
    % Split plant into minimum phase and non-minimum phase parts
    try
        [z, p, k] = zpkdata(G_for_design, 'v');
        
        % Create minimum phase part
        Gp = zpk([], p, k);
        
        % Collect non-minimum phase zeros to be removed
        nmp_zeros = [];
        
        for i = 1:length(z)
            if real(z(i)) <= 0
                % Minimum phase zero - add to Gp
                Gp = zpk([Gp.z; z(i)], Gp.p, Gp.k);
            else
                % Non-minimum phase zero - collect for filter design
                nmp_zeros = [nmp_zeros; z(i)];
            end
        end
        
        if ~isempty(nmp_zeros)
            details = [details, sprintf('Non-minimum phase zeros detected: %d\n', length(nmp_zeros))];
        end
    catch ME
        % If zpk decomposition fails, try a different approach
            warning('designIMC:zpkDecomposition', 'Decomposition failed: %s\nUsing alternative approach.', ME.message);
        
        % Fallback approach using frequency response and approximate FOPDT model
        if ~isnan(plantInfo.FOPDT.K) && ~isnan(plantInfo.FOPDT.T) && ~isnan(plantInfo.FOPDT.L)
            % Use FOPDT model as approximate minimum phase component
            Gp = tf(plantInfo.FOPDT.K, [plantInfo.FOPDT.T, 1]);
            
            % Add delay approximation if significant
            if plantInfo.FOPDT.L > 0.05
                % Pade approximation for delay
                pade_order = 1;
                [num_delay, den_delay] = pade(plantInfo.FOPDT.L, pade_order);
                delay_tf = tf(num_delay, den_delay);
                
                % Combine with main transfer function
                Gp = Gp * delay_tf;
            end
            
            details = [details, 'Using FOPDT model for IMC design due to decomposition failure.\n'];
        else
            % If even FOPDT estimation failed, use plant directly with caution
            Gp = G_for_design;
            details = [details, 'Using original plant for IMC design (non-ideal, may result in instability).\n'];
        end
        
        nmp_zeros = [];
    end
    
    % Adjust lambda based on plant properties for better robustness
    if ~isempty(nmp_zeros)
        % For non-minimum phase plants, increase lambda for robustness
        lambda = lambda * (1 + length(nmp_zeros));
        details = [details, sprintf('Increased lambda to %.4f for non-minimum phase plant.\n', lambda)];
    end
    
    if plantInfo.isHighOrder
        % For high-order systems, increase lambda slightly
        lambda = lambda * 1.2;
        details = [details, sprintf('Increased lambda to %.4f for high-order plant.\n', lambda)];
    end
    
    % IMC filter design based on controller structure
    switch structure
        case 'P'
            % Simple first-order IMC filter
            F = tf(1, [lambda, 1]);
            
        case 'PI'
            % IMC filter for PI-like characteristics
            if plantInfo.hasIntegrator
                % Plant already has integrator - use first-order filter
                F = tf(1, [lambda, 1]);
            else
                % Add integrator to filter
                F = tf(1, [lambda, 1]) * tf(1, [1, 0]);
            end
            
        case 'PD'
            % IMC filter for PD-like characteristics
            F = tf([lambda/4, 1], [lambda, 1]);
            
        case 'PID'
            % IMC filter for PID-like characteristics
            if plantInfo.hasIntegrator
                % Plant already has integrator
                F = tf([lambda/4, 1], [lambda, 1]);
            else
                % Add integrator to filter
                F = tf([lambda/4, 1], [lambda, 1]) * tf(1, [1, 0]);
            end
            
        otherwise
            error('Unsupported controller structure for IMC method');
    end
    
    % Create IMC controller
    try
        % Standard IMC controller formula: Q = F/Gp
        Q = minreal(F / Gp);
        
        % Convert to classical feedback controller: K = Q/(1-Q*G)
        K_raw = minreal(Q / (1 - Q * G_for_design));
    catch ME
        % If direct inversion fails, try numerical approximation
        warning('designIMC:createController','Direct IMC calculation failed: %s\nUsing numerical approximation.', ME.message);
        
        % Get frequency response of plant
        w = logspace(-3, 3, 100);
        [mag_G, phase_G] = bode(G_for_design, w);
        [mag_F, phase_F] = bode(F, w);
        
        % Create approximate Q with bounded gain at high frequencies
        mag_Q = zeros(size(w));
        phase_Q = zeros(size(w));
        
        for i = 1:length(w)
            if mag_G(i) > 1e-3
                mag_Q(i) = mag_F(i) / mag_G(i);
                phase_Q(i) = phase_F(i) - phase_G(i);
            else
                % Limit high-frequency gain
                mag_Q(i) = mag_F(i) / 1e-3;
                phase_Q(i) = phase_F(i) - phase_G(i);
            end
        end
        
        % Create approximate classical controller
        mag_K = zeros(size(w));
        phase_K = zeros(size(w));
        
        for i = 1:length(w)
            den = abs(1 - mag_Q(i) * mag_G(i) * exp(1j * (phase_Q(i) + phase_G(i)) * pi/180));
            mag_K(i) = mag_Q(i) / den;
            phase_K(i) = phase_Q(i) - phase_G(i) * (1 - mag_Q(i) * mag_G(i));
        end
        
        % Fit rational transfer function using frequency response data
        try
            order = min(4, length(plantInfo.poles));
            K_raw = fitsys(w, mag_K, phase_K, order, order);
            details = [details, 'Created approximate controller from frequency response.\n'];
        catch
            % If fitting fails, fall back to a simple controller
            error('Could not fit transfer function to frequency response.');
        end
    end
    
    % Process the controller to match the desired structure
    try
        % Ensure controller is proper and stable
        [num, den] = tfdata(K_raw, 'v');
        
        % Check if controller is proper
        if length(num) > length(den)
            % Add filter to make controller proper
            filter_tf = tf(1, [lambda/10, 1])^(length(num) - length(den));
            K = K_raw * filter_tf;
            details = [details, sprintf('Added filter to make controller proper: (1/(%.4fs+1))^%d\n', ...
                      lambda/10, length(num) - length(den))];
        else
            K = K_raw;
        end
        
        % Check controller stability
        K_poles = pole(K);
        if any(real(K_poles) > 0)
            details = [details, 'Warning: Controller has unstable poles. Applying stabilization.\n'];
            
            % Force controller stability
            [num_K, den_K] = tfdata(K, 'v');
            p = roots(den_K);
            for i = 1:length(p)
                if real(p(i)) > 0
                    p(i) = -abs(real(p(i))) + imag(p(i))*1i;
                end
            end
            den_stable = poly(p);
            K = tf(num_K, den_stable);
        end
        
        % Simplify the controller through model reduction if it's high order
        if length(pole(K)) > 4
            try
                K_simple = balred(K, 4);
                details = [details, 'Controller simplified through model reduction.\n'];
                K = K_simple;
            catch
                % If simplification fails, keep original
                details = [details, 'Controller simplification failed.\n'];
            end
        end
    catch ME
        % If controller processing fails, try creating a structurally correct controller
        warning('designIMC:tryCreateStructurallyCorrectController','Controller processing failed: %s\nCreating a structurally correct controller.', ME.message);
        
        % Create a controller with the desired structure
        switch structure
            case 'P'
                % Simple P controller with gain estimated from IMC principles
                Kp = 1 / plantInfo.FOPDT.K * (plantInfo.FOPDT.T / (lambda + plantInfo.FOPDT.L));
                K = tf(Kp, 1);
                
            case 'PI'
                % PI controller with parameters from IMC design principles
                Kp = 1 / plantInfo.FOPDT.K * (plantInfo.FOPDT.T / (lambda + plantInfo.FOPDT.L));
                Ti = max(plantInfo.FOPDT.T, 4 * lambda);
                K = tf([Kp, Kp/Ti], [1, 0]);
                
            case 'PD'
                % PD controller
                Kp = 1 / plantInfo.FOPDT.K * (plantInfo.FOPDT.T / (lambda + plantInfo.FOPDT.L));
                Td = min(plantInfo.FOPDT.L, lambda/2);
                K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
                
            case 'PID'
                % PID controller
                Kp = 1 / plantInfo.FOPDT.K * (plantInfo.FOPDT.T / (lambda + plantInfo.FOPDT.L));
                Ti = max(plantInfo.FOPDT.T, 4 * lambda);
                Td = min(plantInfo.FOPDT.L, lambda/2);
                K = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
                
            otherwise
                error('Unsupported controller structure');
        end
        
        details = [details, 'Created controller using simplified IMC design rules.\n'];
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
            cl_poles_adj = pole(T_adj);
            
            if all(real(cl_poles_adj) < 0)
                K = K_adjusted;
                details = [details, 'Controller gain reduced by 50% to achieve stability.\n'];
            else
                % Try more aggressive reduction if needed
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