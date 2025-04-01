% LQG design method
function [K, details] = designLQG(G, structure, bandwidth, robustness, epsilon, plantInfo)
    % Enhanced LQG (Linear-Quadratic-Gaussian) controller design
    % Combines optimal LQR state feedback with Kalman filter state estimation
    % Improved handling of plant characteristics and error recovery
    
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    details = [details, sprintf('LQG Design\n----------\nBandwidth target: %.2f rad/s\nRobustness level: %s\n', bandwidth, robustness)];
    
    try
        % Handle unstable plants with special preprocessing
        if plantInfo.isUnstable
            details = [details, 'Plant is unstable. Using modified LQG design approach.\n'];
        end
        
        % Check if Control System Toolbox is available for optimal design
        has_control_toolbox = true;
        if ~exist('lqr', 'file') || ~exist('lqe', 'file')
            has_control_toolbox = false;
            details = [details, 'Control System Toolbox not available. Using approximate LQG design.\n'];
        end
        
        % Convert to state-space for LQG design
        [A, B, C, D] = ssdata(G);
        
        % Basic system analysis
        n = size(A, 1);  % Number of states
        m = size(B, 2);  % Number of inputs
        p = size(C, 1);  % Number of outputs
        
        % Check controllability and observability with robust numerical methods
        try
            % Use SVD-based methods for more reliable analysis
            sv_ctrb = svd(ctrb(A, B));
            sv_obsv = svd(obsv(A, C));
            
            % Check condition of controllability and observability matrices
            cond_ctrb = sv_ctrb(1) / max(sv_ctrb(end), eps);
            cond_obsv = sv_obsv(1) / max(sv_obsv(end), eps);
            
            if cond_ctrb > 1e6
                details = [details, sprintf('Warning: System is poorly controllable (cond = %.1e).\n', cond_ctrb)];
            end
            
            if cond_obsv > 1e6
                details = [details, sprintf('Warning: System is poorly observable (cond = %.1e).\n', cond_obsv)];
            end
            
            % Check rank using tolerance-based approach
            rank_tol = 1e-8 * max(sv_ctrb(1), sv_obsv(1));
            
            is_controllable = (sum(sv_ctrb > rank_tol) == n);
            is_observable = (sum(sv_obsv > rank_tol) == n);
            
            if ~is_controllable
                details = [details, 'Warning: System is not fully controllable. Modifying design approach.\n'];
            end
            
            if ~is_observable
                details = [details, 'Warning: System is not fully observable. Modifying design approach.\n'];
            end
        catch
            % If analysis fails, continue with standard approach
            details = [details, 'Could not perform detailed controllability/observability analysis.\n'];
            is_controllable = true;
            is_observable = true;
        end
        
        % Enhanced weight selection logic based on plant characteristics and requirements
        % Set weighting matrices based on robustness, bandwidth and plant properties
        
        % Base parameters on robustness setting
        switch robustness
            case 'Low'
                % Lower robustness: emphasize performance
                q_scale = 10;    % Higher state penalty
                r_scale = 0.1;   % Lower control penalty
                qn_scale = 10;   % Higher process noise (more aggressive control)
                rn_scale = 1;    % Moderate measurement noise
            case 'High'
                % Higher robustness: more conservative control
                q_scale = 1;     % Moderate state penalty
                r_scale = 10;    % Higher control penalty
                qn_scale = 1;    % Moderate process noise
                rn_scale = 10;   % Higher measurement noise (more filtering)
            otherwise
                % Medium robustness: balanced approach
                q_scale = 5;     % Balanced state penalty
                r_scale = 1;     % Balanced control penalty
                qn_scale = 5;    % Moderate process noise
                rn_scale = 1;    % Moderate measurement noise
        end
        
        % Additional adjustments based on plant characteristics
        if plantInfo.hasRHPZeros
            % For non-minimum phase systems, be more conservative
            r_scale = r_scale * 2;    % Increase control penalty
            details = [details, 'Increased control penalty due to non-minimum phase behavior.\n'];
        end
        
        if plantInfo.isUnstable
            % For unstable systems, be more aggressive with stabilization
            q_scale = q_scale * 2;    % Increase state penalty
            details = [details, 'Increased state penalty for unstable system.\n'];
        end
        
        if plantInfo.isHighOrder
            % For high-order systems, add more filtering
            rn_scale = rn_scale * 1.5;  % More filtering for high-order plants
            details = [details, 'Increased measurement noise for high-order system.\n'];
        end
        
        % Construct the state and control weighting matrices
        Q = C' * C * q_scale;            % Penalize outputs
        R = eye(m) * r_scale;            % Control penalty
        
        % Adjust Q to target the desired bandwidth
        Q = Q * (bandwidth^2);
        
        % Noise covariance matrices
        Qn = eye(n) * qn_scale;          % Process noise
        Rn = eye(p) * rn_scale;          % Measurement noise
        
        % LQG Design approach depends on system characteristics
        if has_control_toolbox && is_controllable && is_observable
            % Standard LQG design if system is well-behaved
            % Design LQR controller (state feedback)
            try
                [K_lqr, S, e] = lqr(A, B, Q, R);
                
                % Design Kalman filter (state estimator)
                [Kf, P, E] = lqe(A, eye(n), C, Qn, Rn);
                
                % Compute LQG controller (combines LQR and Kalman filter)
                Ac = A - B*K_lqr - Kf*C;
                Bc = Kf;
                Cc = -K_lqr;
                Dc = 0;
                
                K_lqg = ss(Ac, Bc, Cc, Dc);
                
                details = [details, sprintf('\nStandard LQG design successful.\n')];
                details = [details, sprintf('LQR gain: [%s]\nKalman gain: [%s]\n', ...
                         mat2str(K_lqr, 3), mat2str(Kf', 3))];
            catch ME
                % If standard LQG fails, use regularized design
                details = [details, sprintf('Standard LQG design failed: %s\nUsing regularized approach.\n', ME.message)];
                
                % Add regularization terms to avoid ill-conditioning
                A_reg = A + 1e-6 * eye(n);
                Q_reg = Q + 1e-6 * eye(n);
                R_reg = R + 1e-6 * eye(m);
                
                % Retry with regularized system
                [K_lqr, ~, ~] = lqr(A_reg, B, Q_reg, R_reg);
                [Kf, ~, ~] = lqe(A_reg, eye(n), C, Qn + 1e-6 * eye(n), Rn + 1e-6 * eye(p));
                
                Ac = A - B*K_lqr - Kf*C;
                Bc = Kf;
                Cc = -K_lqr;
                Dc = 0;
                
                K_lqg = ss(Ac, Bc, Cc, Dc);
                details = [details, 'Regularized LQG design successful.\n'];
            end
        else
            % Alternative approach for problematic systems or missing toolbox
            details = [details, 'Using alternative design approach for problematic system.\n'];
            
            % Create a frequency-domain approximation to LQG
            w = logspace(-3, 3, 100);
            
            try
                % Try loop-shaping approach to approximate LQG behavior
                % Calculate approximate crossover frequency based on bandwidth
                w_c = bandwidth;
                
                % Create a PI filter with phase lead for basic LQG behavior
                Kp = 1 / abs(evalfr(G, 1j*w_c));  % Gain for unity crossover
                Ti = 5 / w_c;  % Integral time constant
                Td = 0.1 / w_c;  % Derivative time constant
                
                K_lqg = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
                
                details = [details, 'Created approximate controller using loop-shaping.\n'];
            catch
                % If even approximate design fails, create a simple controller
                dc_gain = dcgain(G);
                if isnan(dc_gain) || dc_gain == 0 || isinf(dc_gain)
                    dc_gain = 1;  % Default if DC gain is problematic
                end
                
                Kp = 1 / abs(dc_gain);
                Ti = 10;  % Conservative integral action
                Td = 0.1;  % Conservative derivative action
                
                K_lqg = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
                details = [details, 'Created simplified controller due to design difficulties.\n'];
            end
        end
        
        % Convert to desired controller structure with advanced extraction method
        details = [details, sprintf('\nConverting to %s controller structure...\n', structure)];
        
        switch structure
            case 'P'
                % P controller approximation using balanced frequency response
                try
                    % Use mid-frequency gain for better performance
                    w_mid = bandwidth;
                    mag_mid = abs(evalfr(K_lqg, 1j*w_mid));
                    Kp = mag_mid;
                catch
                    % Fallback to DC gain
                    Kp = abs(dcgain(K_lqg));
                    if isnan(Kp) || isinf(Kp)
                        Kp = 1.0;  % Default if computation fails
                    end
                end
                
                K = tf(Kp, 1);
                details = [details, sprintf('P controller with gain: Kp = %.4f', Kp)];
                
            case 'PI'
                % Enhanced PI approximation from LQG controller
                try
                    % Frequency domain approximation for PI controllers
                    w = logspace(-3, 2, 100);
                    [mag, phase] = bode(K_lqg, w);
                    mag = squeeze(mag);
                    phase = squeeze(phase);
                    
                    % Look for integrator (phase near -90° at low frequencies)
                    has_integrator = any(phase(1:min(5, length(phase))) < -70);
                    
                    if has_integrator
                        % Extract parameters from frequency response
                        idx_mid = floor(length(w)/2);
                        Kp = mag(idx_mid);  % Gain near crossover
                        Ki = w(1) * mag(1);  % Integral action from low frequencies
                    else
                        % If no integrator detected, use standard approach
                        Kp = abs(dcgain(K_lqg));
                        if isnan(Kp) || isinf(Kp)
                            Kp = 1.0;
                        end
                        Ki = Kp / 10;  % Conservative integral action
                    end
                    
                    % Adjust for problematic plants
                    if plantInfo.hasIntegrator
                        Ki = Ki * 0.5;  % Reduce integral action for plants with integrator
                    end
                    
                    if plantInfo.isUnstable && plantInfo.stabilizingGain > 0
                        Kp = max(Kp, plantInfo.stabilizingGain * 1.2);  % Ensure stability
                    end
                catch
                    % Fallback to simple design
                    Kp = 1.0;
                    Ki = 0.5;
                end
                
                K = tf([Kp, Ki], [1, 0]);
                details = [details, sprintf('PI controller with:\nKp = %.4f\nKi = %.4f\nTi = %.4f', Kp, Ki, Kp/Ki)];
                
            case 'PD'
                % Enhanced PD extraction with better high-frequency handling
                try
                    % Use frequency response to extract PD parameters
                    w = logspace(-2, 3, 100);
                    [mag, phase] = bode(K_lqg, w);
                    mag = squeeze(mag);
                    phase = squeeze(phase);
                    
                    % Look for derivative action (phase lead at high frequencies)
                    if any(phase > 10)
                        % Find region with maximum phase lead
                        [~, idx_max] = max(phase);
                        w_max = w(idx_max);
                        
                        % Calculate time constants
                        Td = 1 / w_max;
                        
                        % Gain at mid-frequencies
                        idx_mid = floor(length(w)/2);
                        Kp = mag(idx_mid);
                    else
                        % If no clear phase lead, use conservative estimation
                        Kp = abs(dcgain(K_lqg));
                        if isnan(Kp) || isinf(Kp)
                            Kp = 1.0;
                        end
                        Td = 0.1 / bandwidth;  % Based on bandwidth
                    end
                    
                    Kd = Kp * Td;
                    
                    % Adjust for problematic plants
                    if plantInfo.hasRHPZeros
                        Kd = Kd * 0.7;  % Reduce derivative action for non-minimum phase
                    end
                    
                    if plantInfo.isUnstable && plantInfo.stabilizingGain > 0
                        Kp = max(Kp, plantInfo.stabilizingGain * 1.2);  % Ensure stability
                    end
                catch
                    % Fallback to conservative parameters
                    Kp = 1.0;
                    Kd = 0.1;
                end
                
                K = tf([Kd, Kp], [epsilon*Td, 1]);
                details = [details, sprintf('PD controller with:\nKp = %.4f\nKd = %.4f\nTd = %.4f\nFilter epsilon = %.4f', ...
                         Kp, Kd, Kd/Kp, epsilon)];
                
            case 'PID'
                % Advanced PID extraction from full LQG controller
                try
                    % Use frequency response to extract PID parameters
                    w = logspace(-3, 3, 100);
                    [mag, phase] = bode(K_lqg, w);
                    mag = squeeze(mag);
                    phase = squeeze(phase);
                    
                    % Look for both integrator and derivative action
                    has_integrator = any(phase(1:min(5, length(phase))) < -70);
                    has_derivative = any(phase(max(1, end-10):end) > 10);
                    
                    if has_integrator && has_derivative
                        % Full PID behavior detected
                        details = [details, 'PID behavior detected in LQG controller.\n'];
                        
                        % Proportional gain at mid-frequencies
                        idx_mid = floor(length(w)/2);
                        Kp = mag(idx_mid);
                        
                        % Integral gain from low-frequency behavior
                        Ki = w(1) * mag(1);  % Extract from low-frequency gain
                        
                        % Derivative gain from high-frequency behavior
                        [max_phase, idx_max] = max(phase);
                        w_max = w(idx_max);
                        Td = 1 / w_max;
                        Kd = Kp * Td;
                    else
                        % Partial or no PID behavior - conservative approximation
                        details = [details, 'Incomplete PID behavior detected. Using approximation.\n'];
                        
                        % Use DC gain for proportional part
                        Kp = abs(dcgain(K_lqg));
                        if isnan(Kp) || isinf(Kp)
                            Kp = 1.0;
                        end
                        
                        % Conservative time constants based on bandwidth
                        Ti = 10 / bandwidth;
                        Td = 0.1 / bandwidth;
                        
                        Ki = Kp / Ti;
                        Kd = Kp * Td;
                    end
                    
                    % Adjust parameters for problematic plants
                    if plantInfo.hasIntegrator
                        Ki = Ki * 0.5;  % Reduce integral action for plants with integrator
                    end
                    
                    if plantInfo.hasRHPZeros
                        Kd = Kd * 0.7;  % Reduce derivative action for non-minimum phase
                    end
                    
                    if plantInfo.isUnstable && plantInfo.stabilizingGain > 0
                        Kp = max(Kp, plantInfo.stabilizingGain * 1.2);  % Ensure stability
                    end
                catch
                    % Fallback to conservative parameters
                    Kp = 1.0;
                    Ki = 0.1;
                    Kd = 0.1;
                end
                
                % Create filtered PID controller
                K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                
                details = [details, sprintf('PID controller with:\nKp = %.4f\nKi = %.4f\nKd = %.4f\nTi = %.4f\nTd = %.4f\nFilter epsilon = %.4f', ...
                         Kp, Ki, Kd, Kp/Ki, Kd/Kp, epsilon)];
                
            otherwise
                error('Unsupported controller structure for LQG method');
        end
        
        % Verify controller stability
        try
            K_poles = pole(K);
            if any(real(K_poles) > 0)
                details = [details, '\nWarning: Controller has unstable poles. Applying stabilization.\n'];
                
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
        
        % Verify closed-loop stability and adjust if needed
        try
            T = feedback(G*K, 1);
            cl_poles = pole(T);
            
            if any(real(cl_poles) > 0)
                details = [details, '\nWarning: Closed-loop system is unstable. Adjusting controller.\n'];
                
                % Try to stabilize by reducing gain
                [num, den] = tfdata(K, 'v');
                K_adjusted = tf(num * 0.5, den);
                
                % Check if modification helps
                T_adj = feedback(G * K_adjusted, 1);
                
                if all(real(pole(T_adj)) < 0)
                    K = K_adjusted;
                    details = [details, 'Controller gain reduced by 50% to achieve stability.\n'];
                else
                    % Try more aggressive reduction
                    K_adjusted = tf(num * 0.2, den);
                    T_adj = feedback(G * K_adjusted, 1);
                    
                    if all(real(pole(T_adj)) < 0)
                        K = K_adjusted;
                        details = [details, 'Controller gain reduced by 80% to achieve stability.\n'];
                    else
                        details = [details, 'Could not stabilize system by gain reduction. Consider a different method.\n'];
                    end
                end
            else
                % Performance metrics for stable system
                try
                    [Gm, Pm] = margin(G*K);
                    details = [details, sprintf('\nClosed-loop performance:\nGain Margin: %.2f dB\nPhase Margin: %.2f degrees', ...
                             20*log10(Gm), Pm)];
                    
                    % Time-domain performance analysis
                    try
                        info = stepinfo(T);
                        details = [details, sprintf('\nSettling Time: %.2f s\nOvershoot: %.2f%%', ...
                                 info.SettlingTime, info.Overshoot)];
                        
                        % If overshoot is excessive, reduce controller gain
                        if info.Overshoot > 40
                            [num, den] = tfdata(K, 'v');
                            K_adjusted = tf(num * 0.7, den);
                            T_adj = feedback(G * K_adjusted, 1);
                            
                            info_adj = stepinfo(T_adj);
                            if info_adj.Overshoot < info.Overshoot
                                K = K_adjusted;
                                details = [details, sprintf('\nReduced controller gain by 30%% to limit overshoot from %.1f%% to %.1f%%', ...
                                         info.Overshoot, info_adj.Overshoot)];
                            end
                        end
                    catch
                        details = [details, '\nCould not calculate time-domain performance metrics.'];
                    end
                catch
                    details = [details, '\nCould not compute stability margins.'];
                end
            end
        catch
            details = [details, '\nCould not verify closed-loop stability.'];
        end
        
    catch ME
        % Enhanced fallback approach based on plant characteristics
        warning('LQG design error: %s', ME.message);
        details = [details, sprintf('\nLQG design failed: %s\n', ME.message)];
        
        % Create adaptive fallback controller based on plant analysis
        details = [details, 'Using adaptive fallback controller based on plant analysis.\n'];
        
        % Determine appropriate fallback strategy
        if plantInfo.isUnstable
            % For unstable plants, use stabilizing controller with structure-specific modifications
            switch structure
                case 'P'
                    Kp = plantInfo.stabilizingGain * 1.5;
                    K = tf(Kp, 1);
                    details = [details, sprintf('Stabilizing P controller with Kp = %.4f', Kp)];
                    
                case 'PI'
                    Kp = plantInfo.stabilizingGain * 1.5;
                    Ki = Kp * 0.05;  % Very small integral action for stability
                    K = tf([Kp, Ki], [1, 0]);
                    details = [details, sprintf('Stabilizing PI controller with Kp = %.4f, Ki = %.4f', Kp, Ki)];
                    
                case 'PD'
                    Kp = plantInfo.stabilizingGain * 1.2;
                    Td = 0.1;
                    Kd = Kp * Td;
                    K = tf([Kd, Kp], [epsilon*Td, 1]);
                    details = [details, sprintf('Stabilizing PD controller with Kp = %.4f, Kd = %.4f', Kp, Kd)];
                    
                case 'PID'
                    Kp = plantInfo.stabilizingGain * 1.2;
                    Ki = Kp * 0.05;  % Small integral action
                    Td = 0.1;
                    Kd = Kp * Td;
                    K = tf([Kd, Kp, Ki], [epsilon*Td, 1, 0]);
                    details = [details, sprintf('Stabilizing PID controller with Kp = %.4f, Ki = %.4f, Kd = %.4f', Kp, Ki, Kd)];
                    
                otherwise
                    K = tf(plantInfo.stabilizingGain * 1.5, 1);
                    details = [details, 'Using generic stabilizing controller.'];
            end
        else
            % For stable plants, use conservative tuning based on FOPDT approximation
            if ~isnan(plantInfo.FOPDT.K) && ~isnan(plantInfo.FOPDT.T) && ~isnan(plantInfo.FOPDT.L)
                % Use FOPDT parameters for controller design
                Ks = plantInfo.FOPDT.K;
                T = plantInfo.FOPDT.T;
                L = plantInfo.FOPDT.L;
                
                switch structure
                    case 'P'
                        Kp = 0.5 * T / (Ks * L);
                        K = tf(Kp, 1);
                        details = [details, sprintf('Conservative P controller with Kp = %.4f', Kp)];
                        
                    case 'PI'
                        Kp = 0.4 * T / (Ks * L);
                        Ti = 1.2 * T;
                        Ki = Kp / Ti;
                        K = tf([Kp, Ki], [1, 0]);
                        details = [details, sprintf('Conservative PI controller with Kp = %.4f, Ki = %.4f', Kp, Ki)];
                        
                    case 'PD'
                        Kp = 0.5 * T / (Ks * L);
                        Td = 0.5 * L;
                        Kd = Kp * Td;
                        K = tf([Kd, Kp], [epsilon*Td, 1]);
                        details = [details, sprintf('Conservative PD controller with Kp = %.4f, Kd = %.4f', Kp, Kd)];
                        
                    case 'PID'
                        Kp = 0.6 * T / (Ks * L);
                        Ti = 1.0 * T;
                        Td = 0.5 * L;
                        Ki = Kp / Ti;
                        Kd = Kp * Td;
                        K = tf([Kd, Kp, Ki], [epsilon*Td, 1, 0]);
                        details = [details, sprintf('Conservative PID controller with Kp = %.4f, Ki = %.4f, Kd = %.4f', Kp, Ki, Kd)];
                        
                    otherwise
                        K = tf(0.2 / Ks, 1);
                        details = [details, 'Using generic conservative controller.'];
                end
            else
                % Most conservative approach when even FOPDT estimation failed
                dc_gain = abs(plantInfo.dcGain);
                if isnan(dc_gain) || dc_gain == 0 || isinf(dc_gain)
                    dc_gain = 1;
                end
                
                switch structure
                    case 'P'
                        Kp = 0.2 / dc_gain;
                        K = tf(Kp, 1);
                    case 'PI'
                        Kp = 0.2 / dc_gain;
                        Ki = Kp * 0.1;
                        K = tf([Kp, Ki], [1, 0]);
                    case 'PD'
                        Kp = 0.2 / dc_gain;
                        Kd = Kp * 0.1;
                        K = tf([Kd, Kp], [epsilon*Td, 1]);
                    case 'PID'
                        Kp = 0.2 / dc_gain;
                        Ki = Kp * 0.1;
                        Kd = Kp * 0.1;
                        K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                    otherwise
                        K = tf(0.2 / dc_gain, 1);
                end
                
                details = [details, 'Using highly conservative controller due to limited plant information.'];
            end
        end
    end
end