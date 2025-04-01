function [K, details] = designLoopShaping(G, structure, phaseMargin, bandwidth, epsilon, plantInfo)
    % Enhanced Loop-Shaping method with better handling of unstable and difficult plants
    
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    details = [details, sprintf('Desired Phase Margin: %.2f°\nTarget Bandwidth: %.2f rad/s\n', phaseMargin, bandwidth)];
    
    % For unstable plants, we can directly apply loop-shaping
    originalG = G;
    
    % Calculate Bode diagram of the plant
    w = logspace(-3, 4, 1000);
    [mag, phase, wout] = bode(G, w);
    mag = squeeze(mag);
    phase = squeeze(phase);
    
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
            if plantInfo.isUnstable
                % For unstable plants, may need higher gain to reach bandwidth
                Kp = 1.5/mag_at_bw;
            else
                Kp = 1/mag_at_bw;  % Gain for 0dB at bandwidth
            end
            K = tf(Kp, 1);
            
        case 'PI'
            % PI-controller: Kp(1 + 1/(Ti*s))
            Kp = 1/mag_at_bw * 0.8;  % Slightly reduced due to PI phase lag
            
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
            
            K = tf([Kp, Kp/Ti], [1, 0]);
            
        case 'PD'
            % PD-controller: Kp(1 + Td*s)
            if required_phase > 0
                % Calculate Td for desired phase boost
                alpha = (1 - sin(required_phase*pi/180)) / (1 + sin(required_phase*pi/180));
                Td = 1/(bandwidth * sqrt(alpha));
                
                % Calculate Kp for 0dB at bandwidth with phase boost
                [mag_pd, ~] = bode(tf([Td, 1], [epsilon*Td, 1]), bandwidth);
                mag_pd = squeeze(mag_pd);
                Kp = 1 / (mag_at_bw * mag_pd);
            else
                % No phase boost needed
                Kp = 1/mag_at_bw;
                Td = 0.1/bandwidth;
            end
            
            % For unstable plants, increase gain slightly
            if plantInfo.isUnstable
                Kp = Kp * 1.2;
                details = [details, 'Increased gain for unstable plant.\n'];
            end
            
            % For plants with RHP zeros, reduce derivative action
            if plantInfo.hasRHPZeros
                Td = Td * 0.7;
                details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
            end
            
            K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
            
        case 'PID'
            % PID-controller: Kp(1 + 1/(Ti*s) + Td*s)
            if required_phase > 0
                % Calculate Td for phase boost (adjusted for PID)
                alpha = (1 - sin(required_phase*pi/180)) / (1 + sin(required_phase*pi/180));
                Td = 1/(bandwidth * sqrt(alpha));
                
                % Calculate Kp for 0dB at bandwidth with phase boost
                [mag_pd, ~] = bode(tf([Td, 1], [epsilon*Td, 1]), bandwidth);
                mag_pd = squeeze(mag_pd);
                Kp = 1 / (mag_at_bw * mag_pd) * 0.8;  % Reduced for stability
                
                % Ti depends on plant characteristics
                if plantInfo.hasIntegrator
                    Ti = 10/bandwidth;  % Larger Ti for plants with integrator
                else
                    Ti = 5/bandwidth;   % Ti typically below bandwidth frequency
                end
            else
                % No phase boost needed
                Kp = 1/mag_at_bw * 0.8;
                Td = 0.1/bandwidth;
                Ti = 5/bandwidth;
            end
            
            % For unstable plants, increase gain but add more filtering
            if plantInfo.isUnstable
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
            
            K = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
            
        otherwise
            error('Unsupported controller structure for Loop-Shaping method');
    end
    
    % Verify controller
    try
        L = originalG * K;
        [Gm, Pm, ~, ~] = margin(L);
        pm_achieved = Pm;
        gm_achieved = 20*log10(Gm);
        
        details = [details, sprintf('Achieved Phase Margin: %.2f° (Target: %.2f°)\nGain Margin: %.2f dB\n', ...
                  pm_achieved, phaseMargin, gm_achieved)];
        
        % Adjust controller if margins are insufficient
        if pm_achieved < 0.7 * phaseMargin || gm_achieved < 6
            details = [details, 'Insufficient stability margins. Adjusting controller.\n'];
            
            % Reduce gain for better stability
            [num, den] = tfdata(K, 'v');
            K_adjusted = tf(num * 0.7, den);
            
            % Verify improvement
            L_adjusted = originalG * K_adjusted;
            [Gm_adj, Pm_adj] = margin(L_adjusted);
            
            if Pm_adj > pm_achieved && 20*log10(Gm_adj) > gm_achieved
                K = K_adjusted;
                details = [details, sprintf('Controller gain reduced by 30%% for better stability.\n')];
                details = [details, sprintf('New Phase Margin: %.2f°\nNew Gain Margin: %.2f dB\n', ...
                          Pm_adj, 20*log10(Gm_adj))];
            end
        end
    catch
        details = [details, 'Could not compute stability margins.\n'];
    end
    
    % Verify closed-loop stability
    try
        T = feedback(originalG * K, 1);
        if any(real(pole(T)) > 0)
            details = [details, 'Warning: Closed-loop system is unstable. Adding stabilization.\n'];
            
            % Try to stabilize by adjusting controller
            [num, den] = tfdata(K, 'v');
            
            % Reduce gain progressively until stable
            for factor = [0.5, 0.3, 0.1]
                K_test = tf(num * factor, den);
                T_test = feedback(originalG * K_test, 1);
                
                if all(real(pole(T_test)) < 0)
                    K = K_test;
                    details = [details, sprintf('Reduced gain by %.0f%% to achieve stability.\n', (1-factor)*100)];
                    break;
                end
            end
        end
    catch
        details = [details, 'Could not verify closed-loop stability.\n'];
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