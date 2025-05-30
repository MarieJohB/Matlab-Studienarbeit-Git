function [K, details] = designZieglerNicholsOscillation(G, structure, options, plantInfo)
  % Enhanced Ziegler-Nichols Oscillation method with better handling of unstable systems
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % Extract needed parameters
    epsilon = options.epsilon;
    
    % For unstable systems, first pre-stabilize
    if plantInfo.isUnstable
        details = [details, 'Plant is unstable. Using pre-stabilization before Ziegler-Nichols method.\n'];
        
        % Create a robust stabilizing controller for the unstable plant
        K_stab = designRobustStabilizingController(G, plantInfo);
        
        % Create a stabilized version of the plant
        G_stab = feedback(G * K_stab, 1);
        
        % Now perform Ziegler-Nichols analysis on the stabilized plant
        G_for_design = G_stab;
        isPreStabilized = true;
    else
        G_for_design = G;
        isPreStabilized = false;
    end
    
    % Apply the traditional Ziegler-Nichols oscillation method
    try
        % Calculate critical gain and frequency
        Ku = 0;
        Tu = 0;
        
        % Use frequency response method to find critical gain
        w = logspace(-3, 3, 1000);
        [mag, phase] = bode(G_for_design, w);
        mag = squeeze(mag);
        phase = squeeze(phase);
        
        % Find where phase is -180 degrees
        phase_cross_idx = find(phase <= -180, 1);
        
        if ~isempty(phase_cross_idx)
            % Find critical gain and period
            Ku = 1 / mag(phase_cross_idx);
            critical_freq = w(phase_cross_idx);
            Tu = 2 * pi / critical_freq;
            
            details = [details, sprintf('Critical gain Ku = %.4f\n', Ku)];
            details = [details, sprintf('Critical period Tu = %.4f s\n', Tu)];
        else
            % If phase never crosses -180, use approximate method
            [min_phase_val, min_phase_idx] = min(phase);
            phase_margin = 180 + min_phase_val;
            
            if phase_margin < 90
                % Use phase margin to estimate critical gain
                Ku = 1 / mag(min_phase_idx) * (1 + sind(phase_margin)) / sind(phase_margin);
                critical_freq = w(min_phase_idx);
                Tu = 2 * pi / critical_freq;
                
                details = [details, sprintf('Estimated critical gain Ku = %.4f (from phase margin)\n', Ku)];
                details = [details, sprintf('Estimated critical period Tu = %.4f s\n', Tu)];
            else
                error('System does not have sufficient phase lag for Ziegler-Nichols method');
            end
        end
        
        % Calculate controller parameters using ZN rules
        switch structure
            case 'P'
                Kp = 0.5 * Ku;
                K_zn = tf(Kp, 1);
                details = [details, sprintf('P controller: Kp = %.4f\n', Kp)];
                
            case 'PI'
                Kp = 0.45 * Ku;
                Ti = 0.85 * Tu;
                Ki = Kp / Ti;
                K_zn = tf([Kp, Ki], [1, 0]);
                details = [details, sprintf('PI controller: Kp = %.4f, Ti = %.4f, Ki = %.4f\n', Kp, Ti, Ki)];
                
            case 'PD'
                Kp = 0.8 * Ku;
                Td = Tu / 8;
                Kd = Kp * Td;
                K_zn = tf([Kd, Kp], [epsilon*Td, 1]);
                details = [details, sprintf('PD controller: Kp = %.4f, Td = %.4f, Kd = %.4f\n', Kp, Td, Kd)];
                
            case 'PID'
                Kp = 0.6 * Ku;
                Ti = 0.5 * Tu;
                Td = 0.125 * Tu;
                Ki = Kp / Ti;
                Kd = Kp * Td;
                K_zn = tf([Kd, Kp, Ki], [epsilon*Td, 1, 0]);
                details = [details, sprintf('PID controller: Kp = %.4f, Ti = %.4f, Td = %.4f, Ki = %.4f, Kd = %.4f\n', ...
                          Kp, Ti, Td, Ki, Kd)];
                
            otherwise
                error('Unsupported controller structure for Ziegler-Nichols method');
        end
        
        % Adjust gains for better robustness (ZN tends to be aggressive)
        % Apply correction factors based on plant characteristics
        if plantInfo.hasRHPZeros
            details = [details, 'Non-minimum phase behavior detected. Reducing controller gains for stability.\n'];
            [num, den] = tfdata(K_zn, 'v');
            K_zn = tf(num * 0.7, den);  % Reduce gain by 30%
        end
        
        if plantInfo.isHighOrder
            details = [details, 'High-order system detected. Reducing controller gains for stability.\n'];
            [num, den] = tfdata(K_zn, 'v');
            K_zn = tf(num * 0.8, den);  % Reduce gain by 20%
        end
        
        % For pre-stabilized systems, combine with the stabilizing controller
        if isPreStabilized
            K_combined = series(K_zn, K_stab);
            
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
            K = K_zn;
        end
        
    catch ME
        details = [details, sprintf('Error in Ziegler-Nichols oscillation method: %s\n', ME.message)];
        
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