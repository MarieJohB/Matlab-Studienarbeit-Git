function [K, details] = designCHR(G, structure, options, plantInfo)
%% Enhanced controller design methods with improved handling for unstable systems
% This file contains enhanced versions of all controller design methods
% with specific improvements for handling unstable high-order systems
% and better integration with state space models from the app.
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

function [K, details] = enhancedZNStep(G, structure, options, plantInfo)
    % Enhanced Ziegler-Nichols Step Response method with pre-stabilization for unstable plants
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % Extract needed parameters
    epsilon = options.epsilon;
    
    % For unstable systems, first pre-stabilize
    if plantInfo.isUnstable
        details = [details, 'Plant is unstable. Using pre-stabilization before Ziegler-Nichols Step method.\n'];
        
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
    
    % Apply the traditional Ziegler-Nichols step response method
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
        
        % Use tangent method to find delay and time constant
        dy = diff(y_norm) ./ diff(t);
        [max_slope, idx_max_slope] = max(dy);
        
        if max_slope <= 0
            error('Could not determine maximum slope');
        end
        
        % Calculate tangent parameters: y = m*t + b
        m = max_slope;
        b = y_norm(idx_max_slope) - m * t(idx_max_slope);
        
        % Tangent intersections
        t_0 = -b / m;  % Intersection with y=0
        t_1 = (1 - b) / m;  % Intersection with y=1
        
        % Dead time L and rise time T
        L = t_0;
        T = t_1 - t_0;
        K_process = y_final;
        
        if L < 0 || T <= 0
            % Fallback to simpler method if tangent method fails
            idx_28 = find(y_norm >= 0.283, 1);
            idx_63 = find(y_norm >= 0.632, 1);
            
            if isempty(idx_28) || isempty(idx_63)
                error('Could not determine time constants from step response');
            end
            
            L = t(idx_28);
            T = t(idx_63) - t(idx_28);
        end
        
        details = [details, sprintf('Process gain K = %.4f\n', K_process)];
        details = [details, sprintf('Dead time L = %.4f s\n', L)];
        details = [details, sprintf('Time constant T = %.4f s\n', T)];
        
        % Calculate controller parameters using ZN step response rules
        switch structure
            case 'P'
                Kp = (0.9 * T) / (K_process * L);
                K_zn = tf(Kp, 1);
                details = [details, sprintf('P controller: Kp = %.4f\n', Kp)];
                
            case 'PI'
                Kp = (0.9 * T) / (K_process * L);
                Ti = 3.33 * L;
                Ki = Kp / Ti;
                K_zn = tf([Kp, Ki], [1, 0]);
                details = [details, sprintf('PI controller: Kp = %.4f, Ti = %.4f, Ki = %.4f\n', Kp, Ti, Ki)];
                
            case 'PD'
                Kp = (1.2 * T) / (K_process * L);
                Td = 0.5 * L;
                Kd = Kp * Td;
                K_zn = tf([Kd, Kp], [epsilon*Td, 1]);
                details = [details, sprintf('PD controller: Kp = %.4f, Td = %.4f, Kd = %.4f\n', Kp, Td, Kd)];
                
            case 'PID'
                Kp = (1.2 * T) / (K_process * L);
                Ti = 2 * L;
                Td = 0.5 * L;
                Ki = Kp / Ti;
                Kd = Kp * Td;
                K_zn = tf([Kd, Kp, Ki], [epsilon*Td, 1, 0]);
                details = [details, sprintf('PID controller: Kp = %.4f, Ti = %.4f, Td = %.4f, Ki = %.4f, Kd = %.4f\n', ...
                          Kp, Ti, Td, Ki, Kd)];
                
            otherwise
                error('Unsupported controller structure for Ziegler-Nichols step method');
        end
        
        % Apply correction factors for better robustness
        if T/L < 3
            details = [details, 'Low T/L ratio detected. Reducing controller gains for stability.\n'];
            [num, den] = tfdata(K_zn, 'v');
            K_zn = tf(num * 0.7, den);  % Reduce gain by 30%
        end
        
        if plantInfo.hasRHPZeros
            details = [details, 'Non-minimum phase behavior detected. Reducing controller gains for stability.\n'];
            [num, den] = tfdata(K_zn, 'v');
            K_zn = tf(num * 0.7, den);  % Reduce gain by 30%
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
        details = [details, sprintf('Error in Ziegler-Nichols step method: %s\n', ME.message)];
        
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

function [K, details] = enhancedAstrom(G, structure, options, plantInfo)
    % Enhanced Astrom method with pre-stabilization for unstable plants
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % Extract needed parameters
    epsilon = options.epsilon;
    
    % For unstable systems, first pre-stabilize
    if plantInfo.isUnstable
        details = [details, 'Plant is unstable. Using pre-stabilization before Astrom method.\n'];
        
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
    
    % Apply the modified Astrom method
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
                error('System does not have sufficient phase lag for Astrom method');
            end
        end
        
        % Calculate controller parameters using Astrom's improved rules
        % (Modified ZN rules with better robustness)
        switch structure
            case 'P'
                Kp = 0.33 * Ku;  % More conservative than ZN
                K_astrom = tf(Kp, 1);
                details = [details, sprintf('P controller: Kp = %.4f\n', Kp)];
                
            case 'PI'
                Kp = 0.33 * Ku;  % More conservative than ZN
                Ti = 0.9 * Tu;  % Slightly modified from ZN
                Ki = Kp / Ti;
                K_astrom = tf([Kp, Ki], [1, 0]);
                details = [details, sprintf('PI controller: Kp = %.4f, Ti = %.4f, Ki = %.4f\n', Kp, Ti, Ki)];
                
            case 'PD'
                Kp = 0.67 * Ku;  % More conservative than ZN
                Td = Tu / 7;     % Modified from ZN
                Kd = Kp * Td;
                K_astrom = tf([Kd, Kp], [epsilon*Td, 1]);
                details = [details, sprintf('PD controller: Kp = %.4f, Td = %.4f, Kd = %.4f\n', Kp, Td, Kd)];
                
            case 'PID'
                Kp = 0.45 * Ku;  % More conservative than ZN
                Ti = 0.85 * Tu;  % Slightly modified from ZN
                Td = 0.16 * Tu;  % Slightly modified from ZN
                Ki = Kp / Ti;
                Kd = Kp * Td;
                K_astrom = tf([Kd, Kp, Ki], [epsilon*Td, 1, 0]);
                details = [details, sprintf('PID controller: Kp = %.4f, Ti = %.4f, Td = %.4f, Ki = %.4f, Kd = %.4f\n', ...
                          Kp, Ti, Td, Ki, Kd)];
                
            otherwise
                error('Unsupported controller structure for Astrom method');
        end
        
        % Adjust for plant characteristics
        if plantInfo.hasRHPZeros
            details = [details, 'Non-minimum phase behavior detected. Reducing controller gains for stability.\n'];
            [num, den] = tfdata(K_astrom, 'v');
            K_astrom = tf(num * 0.8, den);  % Reduce gain by 20%
        end
        
        % For pre-stabilized systems, combine with the stabilizing controller
        if isPreStabilized
            K_combined = series(K_astrom, K_stab);
            
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
            K = K_astrom;
        end
        
    catch ME
        details = [details, sprintf('Error in Astrom method: %s\n', ME.message)];
        
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

function [K, details] = enhancedCHR(G, structure, options, plantInfo)
    % Enhanced CHR (Chien-Hrones-Reswick) method with pre-stabilization for unstable plants
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % Extract needed parameters
    epsilon = options.epsilon;
    
    % For unstable systems, first pre-stabilize
    if plantInfo.isUnstable
        details = [details, 'Plant is unstable. Using pre-stabilization before CHR method.\n'];
        
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
    
    % Apply the CHR method
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
            K_process = plantInfo.FOPDT.K;
        else
            % Use tangent method to find delay and time constant
            try
                dy = diff(y_norm) ./ diff(t);
                [max_slope, idx_max_slope] = max(dy);
                
                if max_slope <= 0
                    error('Could not determine maximum slope');
                end
                
                % Calculate tangent parameters: y = m*t + b
                m = max_slope;
                b = y_norm(idx_max_slope) - m * t(idx_max_slope);
                
                % Tangent intersections
                t_0 = -b / m;  % Intersection with y=0
                t_1 = (1 - b) / m;  % Intersection with y=1
                
                % Dead time L and rise time T
                L = t_0;
                T = t_1 - t_0;
            catch
                % Fallback to simpler method if tangent method fails
                idx_28 = find(y_norm >= 0.283, 1);
                idx_63 = find(y_norm >= 0.632, 1);
                
                if isempty(idx_28) || isempty(idx_63)
                    error('Could not determine time constants from step response');
                end
                
                L = t(idx_28);
                T = t(idx_63) - t(idx_28);
            end
            
            K_process = y_final;
        end
        
        details = [details, sprintf('Process gain K = %.4f\n', K_process)];
        details = [details, sprintf('Dead time L = %.4f s\n', L)];
        details = [details, sprintf('Time constant T = %.4f s\n', T)];
        
        % CHR rules for 0% overshoot (more robust than ZN)
        switch structure
            case 'P'
                Kp = (0.3 * T) / (K_process * L);
                K_chr = tf(Kp, 1);
                details = [details, sprintf('P controller: Kp = %.4f\n', Kp)];
                
            case 'PI'
                Kp = (0.35 * T) / (K_process * L);
                Ti = 1.16 * T;
                Ki = Kp / Ti;
                K_chr = tf([Kp, Ki], [1, 0]);
                details = [details, sprintf('PI controller: Kp = %.4f, Ti = %.4f, Ki = %.4f\n', Kp, Ti, Ki)];
                
            case 'PD'
                Kp = (0.5 * T) / (K_process * L);
                Td = 0.5 * L;
                Kd = Kp * Td;
                K_chr = tf([Kd, Kp], [epsilon*Td, 1]);
                details = [details, sprintf('PD controller: Kp = %.4f, Td = %.4f, Kd = %.4f\n', Kp, Td, Kd)];
                
            case 'PID'
                Kp = (0.6 * T) / (K_process * L);
                Ti = T;
                Td = 0.5 * L;
                Ki = Kp / Ti;
                Kd = Kp * Td;
                K_chr = tf([Kd, Kp, Ki], [epsilon*Td, 1, 0]);
                details = [details, sprintf('PID controller: Kp = %.4f, Ti = %.4f, Td = %.4f, Ki = %.4f, Kd = %.4f\n', ...
                          Kp, Ti, Td, Ki, Kd)];
                
            otherwise
                error('Unsupported controller structure for CHR method');
        end
        
        % For pre-stabilized systems, combine with the stabilizing controller
        if isPreStabilized
            K_combined = series(K_chr, K_stab);
            
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
            K = K_chr;
        end
        
    catch ME
        details = [details, sprintf('Error in CHR method: %s\n', ME.message)];
        
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