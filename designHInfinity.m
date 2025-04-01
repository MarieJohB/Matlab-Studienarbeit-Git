% H-infinity design method
function [K, details] = designHInfinity(G, structure, robustness, epsilon, plantInfo)
    % Enhanced H-infinity controller design focusing on robust stability and performance
    
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % Map robustness setting to gamma value (smaller = more robust but conservative)
    switch robustness
        case 'Low'
            gamma = 3.0;  % Less conservative
        case 'Medium'
            gamma = 2.0;  % Balanced approach
        case 'High'
            gamma = 1.2;  % More robust
        otherwise
            gamma = 2.0;  % Default
    end
    
    details = [details, sprintf('H-infinity Design\n----------------\nRobustness level: %s\nGamma value: %.2f\n', robustness, gamma)];
    
    try
        % For unstable plants, we don't need pre-stabilization since H-infinity directly handles instability
        % However, we'll analyze the plant to determine appropriate weighting functions
        
        % Convert to state-space for H-infinity synthesis
        [A, B, C, D] = ssdata(G);
        nx = size(A, 1);  % Number of states
        
        % Determine plant bandwidth for weight selection
        try
            w_bandwidth = getBandwidth(G, plantInfo);
            details = [details, sprintf('Estimated plant bandwidth: %.4f rad/s\n', w_bandwidth)];
        catch
            % If bandwidth estimation fails, use a default value
            w_bandwidth = 1.0;
            details = [details, 'Could not determine bandwidth. Using default value.\n'];
        end
        
        % Determine crossover frequency target
        if plantInfo.hasRHPZeros
            % For non-minimum phase systems, limit bandwidth
            z = plantInfo.zeros;
            rhp_zeros = z(real(z) > 0);
            min_rhp_zero = min(abs(rhp_zeros));
            w_c = min(w_bandwidth, 0.5 * min_rhp_zero);
            details = [details, sprintf('Limiting bandwidth due to RHP zero at %.4f rad/s\n', min_rhp_zero)];
        else
            w_c = w_bandwidth;
        end
        
        % Define weighting functions based on desired properties and plant characteristics
        % Use inverse of desired sensitivity/complementary sensitivity shapes
        if strcmpi(robustness, 'High')
            % Higher robustness: emphasize uncertainty rejection
            M_s = 1.2;  % Maximum sensitivity peak
            A_t = 0.01; % Maximum high-frequency gain
            
            % Enhanced weighting function selection for high robustness
            w_b = 0.01 * w_c;  % Low-frequency performance bound
            Ws = tf([1/M_s, w_b], [1, w_b/1000]);  % S should be small at low frequencies
            
            w_bc = 10 * w_c;  % Control bandwidth limit
            Wks = tf([1, w_bc/10], [0.001, w_bc]);  % KS weight (control effort)
            
            Wt = tf([1, w_c], [A_t, w_c*10]);  % T should be small at high frequencies
        elseif strcmpi(robustness, 'Low')
            % Lower robustness: emphasize performance
            M_s = 2.0;  % Higher sensitivity peak allowed
            A_t = 0.1;  % Higher complementary sensitivity at high frequencies
            
            w_b = 0.05 * w_c;  % Low-frequency performance bound
            Ws = tf([1/M_s, w_b], [1, w_b/100]);  % S should be small at low frequencies
            
            w_bc = 20 * w_c;  % Higher control bandwidth
            Wks = tf([0.1, w_bc/5], [0.01, w_bc]);  % Less aggressive control effort limitation
            
            Wt = tf([1, w_c/2], [A_t, w_c*5]);  % T should be small at high frequencies
        else
            % Medium robustness: balanced approach
            M_s = 1.5;  % Moderate sensitivity peak
            A_t = 0.05; % Moderate high-frequency gain
            
            w_b = 0.02 * w_c;  % Low-frequency performance bound
            Ws = tf([1/M_s, w_b], [1, w_b/500]);  % S should be small at low frequencies
            
            w_bc = 15 * w_c;  % Moderate control bandwidth
            Wks = tf([0.5, w_bc/8], [0.005, w_bc]);  % Moderate control effort limitation
            
            Wt = tf([1, w_c/1.5], [A_t, w_c*8]);  % T should be small at high frequencies
        end
        
        % Adjust weights for specific plant characteristics
        if plantInfo.isUnstable
            % For unstable plants, increase the sensitivity weight at crossover
            Ws = Ws * tf([1, w_c/2], [1, w_c/20]);
            details = [details, 'Adjusted weights for unstable plant.\n'];
        end
        
        if plantInfo.hasRHPZeros
            % For non-minimum phase systems, reduce the bandwidth of complementary sensitivity
            Wt = Wt * tf([1, min_rhp_zero/3], [1, min_rhp_zero*2]);
            details = [details, 'Adjusted weights for non-minimum phase behavior.\n'];
        end
        
        if plantInfo.hasIntegrator
            % For plants with integrator, reduce sensitivity weight at low frequencies
            Ws = Ws * tf([1, w_b/10], [1, w_b/100]);
            details = [details, 'Adjusted weights for plant with integrator.\n'];
        end
        
        details = [details, sprintf('\nSensitivity weight Ws: %s\nControl effort weight Wks: %s\nRobustness weight Wt: %s\n', ...
                 char(Ws), char(Wks), char(Wt))];
        
        % For simple cases, we'll use loop-shaping as an approximation to H-infinity synthesis
        details = [details, 'Using loop-shaping approximation for H-infinity design.\n'];
        K_temp = loopsyn(G, Ws);
        
        % Extract fixed-structure controller from H-infinity result
        switch structure
            case 'P'
                % P controller: extract proportional gain from H-infinity controller
                [num, den] = tfdata(K_temp, 'v');
                
                % Extract proportional gain as DC gain
                Kp = dcgain(K_temp);
                
                % Ensure the gain is reasonable based on plant characteristics
                if plantInfo.isUnstable
                    % For unstable plants, ensure gain is sufficient for stabilization
                    Kp = max(Kp, 1.2 * plantInfo.stabilizingGain);
                end
                
                K = tf(Kp, 1);
                details = [details, sprintf('\nP controller with gain: Kp = %.4f', Kp)];
                
            case 'PI'
                % Extract approximate PI parameters from H-infinity controller
                [num, den] = tfdata(K_temp, 'v');
                
                % Convert to frequency domain for characteristics extraction
                w = logspace(-3, 3, 100);
                [mag, phase] = bode(K_temp, w);
                mag = squeeze(mag);
                phase = squeeze(phase);
                
                % Estimate PI parameters from frequency response
                % Proportional gain near crossover frequency
                idx_mid = floor(length(w)/2);
                Kp = mag(idx_mid);
                
                % Integral gain from low-frequency behavior
                phase_low = phase(1);
                if phase_low < -85  % Close to -90 degrees indicates integrator
                    Ki = w(1) * mag(1);  % Ki ≈ ω * |K(jω)| for ω→0
                else
                    % Estimate from the rise in gain at low frequencies
                    if length(w) > 10
                        slope_low = (log10(mag(5)) - log10(mag(1))) / (log10(w(5)) - log10(w(1)));
                        if slope_low < -0.5  % Indicates integral action
                            Ki = w(1) * mag(1);
                        else
                            Ki = Kp * 0.1;  % Conservative default
                        end
                    else
                        Ki = Kp * 0.1;  % Conservative default
                    end
                end
                
                % Adjust parameters based on plant characteristics
                if plantInfo.hasIntegrator
                    Ki = Ki * 0.5;  % Reduce integral action for plants with integrator
                    details = [details, 'Reduced integral gain for plant with integrator.\n'];
                end
                
                if plantInfo.isUnstable
                    Kp = max(Kp, 1.2 * plantInfo.stabilizingGain);  % Ensure stability
                    details = [details, 'Adjusted Kp to ensure stabilization.\n'];
                end
                
                K = tf([Kp, Ki], [1, 0]);
                details = [details, sprintf('\nPI controller with:\nKp = %.4f\nKi = %.4f\nTi = %.4f', Kp, Ki, Kp/Ki)];
                
            case 'PD'
                % Extract approximate PD parameters from H-infinity controller
                w = logspace(-3, 3, 100);
                [mag, phase] = bode(K_temp, w);
                mag = squeeze(mag);
                phase = squeeze(phase);
                
                % Look for phase lead (positive phase) indicating derivative action
                if any(phase > 5)
                    % Find region with maximum phase lead
                    [max_phase, idx_max] = max(phase);
                    w_lead = w(idx_max);
                    
                    % Estimate Td based on frequency of maximum phase lead
                    Td = 1 / (2 * w_lead);
                    
                    % Proportional gain at mid-frequencies
                    idx_mid = floor(length(w)/2);
                    Kp = mag(idx_mid);
                    
                    % Derivative gain
                    Kd = Kp * Td;
                else
                    % If no clear phase lead, use conservative parameters
                    Kp = dcgain(K_temp);
                    
                    % For plants with delay, adjust Td based on estimated delay
                    if plantInfo.hasDelay && ~isnan(plantInfo.FOPDT.L)
                        Td = 0.1 * plantInfo.FOPDT.L;
                    else
                        Td = 0.1;  % Default
                    end
                    
                    Kd = Kp * Td;
                    details = [details, 'No clear phase lead detected. Using conservative PD parameters.\n'];
                end
                
                % Adjust parameters for problematic plants
                if plantInfo.hasRHPZeros
                    Td = Td * 0.7;  % Reduce derivative action for non-minimum phase systems
                    details = [details, 'Reduced derivative action for non-minimum phase plant.\n'];
                end
                
                if plantInfo.isUnstable
                    Kp = max(Kp, 1.2 * plantInfo.stabilizingGain);  % Ensure stability
                    details = [details, 'Adjusted Kp to ensure stabilization.\n'];
                end
                
                K = tf([Kd, Kp], [epsilon*Td, 1]);
                details = [details, sprintf('\nPD controller with:\nKp = %.4f\nKd = %.4f\nTd = %.4f\nFilter epsilon = %.4f', ...
                          Kp, Kd, Td, epsilon)];
                
            case 'PID'
                % Extract PID parameters from H-infinity controller frequency response
                w = logspace(-4, 4, 200);
                [mag, phase] = bode(K_temp, w);
                mag = squeeze(mag);
                phase = squeeze(phase);
                
                % Look for integrator (phase approaching -90° at low frequencies)
                phase_low = mean(phase(1:min(5, length(phase))));
                has_integrator = (phase_low < -75);
                
                % Look for derivative action (phase approaching +90° at high frequencies)
                phase_high = mean(phase(max(1, length(phase)-5):end));
                has_derivative = (phase_high > 10);
                
                % Estimate PID parameters from frequency response
                if has_integrator && has_derivative
                    % Full PID behavior detected
                    details = [details, 'Full PID behavior detected in H-infinity controller.\n'];
                    
                    % Find crossover frequency (where phase ≈ 0°)
                    crossover_idx = find(abs(phase) < 30, 1);
                    if isempty(crossover_idx)
                        crossover_idx = floor(length(w)/2);  % Default to middle frequency
                    end
                    w_c = w(crossover_idx);
                    
                    % Proportional gain near crossover
                    Kp = mag(crossover_idx);
                    
                    % Integral gain from low frequency behavior
                    Ki = w(1) * mag(1);  % Ki ≈ ω * |K(jω)| at low frequencies
                    
                    % Derivative gain from high frequency behavior
                    Kd = mag(end) / w(end);  % Approximate for filtered derivative
                else
                    % Partial behavior - estimate gains conservatively
                    details = [details, 'Incomplete PID behavior. Using conservative estimation.\n'];
                    
                    % Get DC gain for proportional term
                    try
                        Kp = dcgain(K_temp);
                    catch
                        Kp = mag(floor(length(mag)/2));  % Use gain at mid-frequency
                    end
                    
                    % Conservative integral and derivative gains
                    Ki = Kp * 0.2;
                    Kd = Kp * 0.1;
                end
                
                % Adjust parameters based on plant characteristics
                if plantInfo.hasIntegrator
                    Ki = Ki * 0.5;  % Reduce integral action
                    details = [details, 'Reduced integral gain for plant with integrator.\n'];
                end
                
                if plantInfo.hasRHPZeros
                    Kd = Kd * 0.7;  % Reduce derivative action
                    details = [details, 'Reduced derivative action for non-minimum phase plant.\n'];
                end
                
                if plantInfo.isUnstable
                    Kp = max(Kp, 1.2 * plantInfo.stabilizingGain);  % Ensure stability
                    details = [details, 'Adjusted Kp to ensure stabilization.\n'];
                end
                
                % Calculate traditional time constants
                Ti = Kp / Ki;
                Td = Kd / Kp;
                
                K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                details = [details, sprintf('\nPID controller with:\nKp = %.4f\nKi = %.4f\nKd = %.4f\nTi = %.4f\nTd = %.4f\nFilter epsilon = %.4f', ...
                          Kp, Ki, Kd, Ti, Td, epsilon)];
                
            otherwise
                error('Unsupported controller structure for H-infinity method');
        end
        
        % Verify controller stability
        try
            K_poles = pole(K);
            if any(real(K_poles) > 0)
                details = [details, '\nWarning: Controller contains unstable poles. Applying stabilization.\n'];
                
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
        
        % Verify closed-loop stability
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
                % Performance metrics if stable
                try
                    [Gm, Pm] = margin(G*K);
                    details = [details, sprintf('\nClosed-loop performance:\nGain Margin: %.2f dB\nPhase Margin: %.2f degrees', ...
                             20*log10(Gm), Pm)];
                catch
                    details = [details, '\nCould not compute stability margins.'];
                end
            end
        catch
            details = [details, '\nCould not verify closed-loop stability.'];
        end
    catch ME
        % Enhanced error handling with more information
        warning('H-infinity design error: %s', ME.message);
        details = [details, sprintf('\nH-infinity design failed: %s\n', ME.message)];
        
        % Create fallback controller based on plant characteristics
        details = [details, 'Using adaptive fallback controller based on plant analysis.\n'];
        
        % Determine appropriate fallback strategy
        if plantInfo.isUnstable
            % For unstable plants, use stabilizing controller
            % with structure-specific modifications
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
                    K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
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
                        K = tf([Kd, Kp], [epsilon*Kd, 1]);
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