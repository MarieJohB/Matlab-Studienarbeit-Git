function [K, details] = designLoopShaping(G, structure, options, plantInfo)
     % Enhanced Loop-Shaping method with improved handling of unstable systems
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % Extract needed parameters
    epsilon = options.epsilon;
    phaseMargin = options.phaseMargin;
    bandwidth = options.bandwidth;
    
    details = [details, sprintf('Desired Phase Margin: %.2f°\nTarget Bandwidth: %.2f rad/s\n', phaseMargin, bandwidth)];
    
    % For unstable plants, consider using state space approach if available
    if plantInfo.isUnstable && isfield(options, 'stateSpace')
        details = [details, 'Unstable plant with state space model available. Using state space approach.\n'];
        
        try
            % Extract state space model
            sys_ss = options.stateSpace;
            
            % Design controller using pole placement with loop-shaping objectives
            opt_pp = options;
            opt_pp.damping = 0.8;  % Good damping for stability
            
            [K, details_pp] = designPolePlacement(G, structure, opt_pp, plantInfo);
            
            details = [details, details_pp];
            return;
        catch ME
            details = [details, sprintf('State space approach failed: %s\n', ME.message)];
            details = [details, 'Falling back to traditional loop-shaping.\n'];
        end
    end
    
    % For unstable plants, we need to pre-stabilize first
    if plantInfo.isUnstable
        details = [details, 'Plant is unstable. Using pre-stabilization before loop-shaping.\n'];
        
        % Create a robust stabilizing controller for the unstable plant
        K_stab = designRobustStabilizingController(G, plantInfo);
        
        % Create a stabilized version of the plant
        G_stab = feedback(G * K_stab, 1);
        G_for_design = G_stab;
        isPreStabilized = true;
    else
        % For stable plants, use direct loop-shaping
        G_for_design = G;
        isPreStabilized = false;
    end
    
    % Detect potential bandwidth limitations
    if plantInfo.hasRHPZeros
        % Find RHP zero frequencies
        z = plantInfo.zeros;
        rhp_zeros = z(real(z) > 0);
        rhp_zero_freqs = abs(rhp_zeros);
        min_rhp_zero_freq = min(rhp_zero_freqs);
        
        % Warn if specified bandwidth exceeds limitations
        if bandwidth > 0.5 * min_rhp_zero_freq
            details = [details, sprintf('Warning: Target bandwidth (%.2f rad/s) exceeds RHP zero limitation (%.2f rad/s).\n', ...
                       bandwidth, 0.5 * min_rhp_zero_freq)];
            % Adjust bandwidth to a safer value
            bandwidth = 0.4 * min_rhp_zero_freq;
            details = [details, sprintf('Adjusting bandwidth to %.2f rad/s for stability.\n', bandwidth)];
        end
    end
    
    % Calculate Bode diagram of the plant
    w = logspace(-3, 4, 1000);
    [mag, phase, wout] = bode(G_for_design, w);
    mag = squeeze(mag);
    phase = squeeze(phase);
    
    % Find index for the desired bandwidth
    [~, idx] = min(abs(wout - bandwidth));
    phase_at_bw = phase(idx);
    mag_at_bw = mag(idx);
    
    % Calculate required phase boost
    required_phase = phaseMargin - (180 + phase_at_bw);
    
    details = [details, sprintf('Phase at bandwidth: %.2f°\nRequired phase boost: %.2f°\n', phase_at_bw, required_phase)];
    
    % Controller design based on structure
    switch structure
        case 'P'
            % Simple P-controller
            if isPreStabilized || plantInfo.isHighOrder
                % For pre-stabilized or high-order plants, use more conservative gain
                Kp = 0.7/mag_at_bw;
            else
                Kp = 1/mag_at_bw;  % Gain for 0dB at bandwidth
            end
            K_ls = tf(Kp, 1);
            details = [details, sprintf('P controller: Kp = %.4f\n', Kp)];
            
        case 'PI'
            % PI-controller: Kp(1 + 1/(Ti*s))
            if isPreStabilized || plantInfo.isHighOrder
                % For pre-stabilized or high-order plants, use more conservative gain
                Kp = 0.7/mag_at_bw;
            else
                Kp = 0.8/mag_at_bw;  % Slightly reduced due to PI phase lag
            end
            
            % Choose Ti based on plant characteristics
            if plantInfo.hasIntegrator
                Ti = 10/bandwidth;  % Larger Ti for plants with integrator
            else
                Ti = 5/bandwidth;   % Ti typically below bandwidth frequency
            end
            
            % For plants with RHP zeros, adjust Ti for stability
            if plantInfo.hasRHPZeros
                Ti = Ti * 1.5;  % Slower integral action
                details = [details, 'Adjusting Ti for non-minimum phase behavior.\n'];
            end
            
            Ki = Kp / Ti;
            K_ls = tf([Kp, Ki], [1, 0]);
            details = [details, sprintf('PI controller: Kp = %.4f, Ti = %.4f, Ki = %.4f\n', Kp, Ti, Ki)];
            
        case 'PD'
            % PD-controller: Kp(1 + Td*s)
            if required_phase > 0
                % Calculate Td for desired phase boost
                alpha = (1 - sin(required_phase*pi/180)) / (1 + sin(required_phase*pi/180));
                Td = 1/(bandwidth * sqrt(alpha));
                
                % Calculate Kp for 0dB at bandwidth with phase boost
                [mag_pd, ~] = bode(tf([Td, 1], [epsilon*Td, 1]), bandwidth);
                mag_pd = squeeze(mag_pd);
                
                if isPreStabilized || plantInfo.isHighOrder
                    Kp = 0.7 / (mag_at_bw * mag_pd);  % More conservative
                else
                    Kp = 1 / (mag_at_bw * mag_pd);
                end
            else
                % No phase boost needed
                if isPreStabilized || plantInfo.isHighOrder
                    Kp = 0.7/mag_at_bw;
                else
                    Kp = 1/mag_at_bw;
                end
                Td = 0.1/bandwidth;
            end
            
            % For unstable plants, increase gain slightly
            if plantInfo.isUnstable && ~isPreStabilized
                Kp = Kp * 1.2;
                details = [details, 'Increased gain for unstable plant.\n'];
            end
            
            % For plants with RHP zeros, reduce derivative action
            if plantInfo.hasRHPZeros
                Td = Td * 0.7;
                details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
            end
            
            Kd = Kp * Td;
            K_ls = tf([Kd, Kp], [epsilon*Td, 1]);
            details = [details, sprintf('PD controller: Kp = %.4f, Td = %.4f, Kd = %.4f\n', Kp, Td, Kd)];
            
        case 'PID'
            % PID-controller: Kp(1 + 1/(Ti*s) + Td*s)
            if required_phase > 0
                % Calculate Td for phase boost (adjusted for PID)
                alpha = (1 - sin(required_phase*pi/180)) / (1 + sin(required_phase*pi/180));
                Td = 1/(bandwidth * sqrt(alpha));
                
                % Calculate Kp for 0dB at bandwidth with phase boost
                [mag_pd, ~] = bode(tf([Td, 1], [epsilon*Td, 1]), bandwidth);
                mag_pd = squeeze(mag_pd);
                
                if isPreStabilized || plantInfo.isHighOrder
                    Kp = 0.6 / (mag_at_bw * mag_pd);  % More conservative
                else
                    Kp = 0.8 / (mag_at_bw * mag_pd);
                end
                
                % Ti depends on plant characteristics
                if plantInfo.hasIntegrator
                    Ti = 10/bandwidth;  % Larger Ti for plants with integrator
                else
                    Ti = 5/bandwidth;   % Ti typically below bandwidth frequency
                end
            else
                % No phase boost needed
                if isPreStabilized || plantInfo.isHighOrder
                    Kp = 0.6/mag_at_bw;
                else
                    Kp = 0.8/mag_at_bw;
                end
                Td = 0.1/bandwidth;
                Ti = 5/bandwidth;
            end
            
            % For unstable plants, increase gain but add more filtering
            if plantInfo.isUnstable && ~isPreStabilized
                Kp = Kp * 1.1;
                epsilon = epsilon * 2;  % More filtering for stability
                details = [details, 'Adjusted parameters for unstable plant.\n'];
            end
            
            % For plants with RHP zeros, reduce derivative action
            if plantInfo.hasRHPZeros
                Td = Td * 0.6;
                Ti = Ti * 1.5;  % Slower integral action
                details = [details, 'Adjusted parameters for non-minimum phase behavior.\n'];
            end
            
            Ki = Kp / Ti;
            Kd = Kp * Td;
            
            K_ls = tf([Kd, Kp, Ki], [epsilon*Td, 1, 0]);
            details = [details, sprintf('PID controller: Kp = %.4f, Ti = %.4f, Td = %.4f, Ki = %.4f, Kd = %.4f\n', ...
                      Kp, Ti, Td, Ki, Kd)];
            
        otherwise
            error('Unsupported controller structure for Loop-Shaping method');
    end
    
    % Verify controller performance
    try
        L = G_for_design * K_ls;
        [Gm, Pm, ~, ~] = margin(L);
        pm_achieved = Pm;
        gm_achieved = 20*log10(Gm);
        
        details = [details, sprintf('Achieved Phase Margin: %.2f° (Target: %.2f°)\nGain Margin: %.2f dB\n', ...
                  pm_achieved, phaseMargin, gm_achieved)];
        
        % Adjust controller if margins are insufficient
        if pm_achieved < 0.7 * phaseMargin || gm_achieved < 6
            details = [details, 'Insufficient stability margins. Adjusting controller.\n'];
            
            % Reduce gain for better stability
            [num, den] = tfdata(K_ls, 'v');
            K_ls = tf(num * 0.7, den);
            
            % Verify improvement
            L_adjusted = G_for_design * K_ls;
            [Gm_adj, Pm_adj] = margin(L_adjusted);
            
            details = [details, sprintf('Controller gain reduced by 30%% for better stability.\n')];
            details = [details, sprintf('New Phase Margin: %.2f°\nNew Gain Margin: %.2f dB\n', ...
                      Pm_adj, 20*log10(Gm_adj))];
        end
    catch ME
        details = [details, sprintf('Could not compute stability margins: %s\n', ME.message)];
    end
    
    % For pre-stabilized systems, combine with the stabilizing controller
    if isPreStabilized
        K_combined = series(K_ls, K_stab);
        
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
        K = K_ls;
    end
    
    % Verify closed-loop stability
    try
        T = feedback(G * K, 1);
        if any(real(pole(T)) > 0)
            details = [details, 'Warning: Closed-loop system is unstable. Adding stabilization.\n'];
            
            % Try to stabilize by adjusting controller
            [num, den] = tfdata(K, 'v');
            
            % Reduce gain progressively until stable
            for factor = [0.5, 0.3, 0.1]
                K_test = tf(num * factor, den);
                T_test = feedback(G * K_test, 1);
                
                if all(real(pole(T_test)) < 0)
                    K = K_test;
                    details = [details, sprintf('Reduced gain by %.0f%% to achieve stability.\n', (1-factor)*100)];
                    break;
                end
            end
        end
    catch ME
        details = [details, sprintf('Could not verify closed-loop stability: %s\n', ME.message)];
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
end