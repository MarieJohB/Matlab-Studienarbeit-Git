function [K, details] = designHInfinity(G, structure, robustness, epsilon, plantInfo)
    % Enhanced H-infinity controller design focusing on robust stability and performance
    % Specifically improved for high-order unstable systems with better state-space integration
    
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % Map robustness setting to gamma value (smaller = more robust but conservative)
    if isstruct(robustness)
        if isfield(robustness, 'robustness')
            robustness = robustness.robustness;
        else
            robustness = 'Medium';
        end
    end
    
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
    
    % Check if input is state-space or transfer function
    isStateSpace = isa(G, 'ss');
    
    % Try to use state-space model if available
    if ~isStateSpace
        try
            if exist('stateSpace', 'var')
                G_ss = stateSpace;
                isStateSpace = true;
                details = [details, 'Using state-space model from options.\n'];
            else
                G_ss = ss(G);
                isStateSpace = true;
                details = [details, 'Successfully converted to state-space representation.\n'];
            end
        catch ME
            details = [details, sprintf('Could not convert to state-space: %s\n', ME.message)];
            G_ss = [];
        end
    else
        G_ss = G;
        details = [details, 'Using provided state-space model directly.\n'];
    end
    
    try
        % For unstable plants, handle differently depending on model type
        if plantInfo.isUnstable
            if isStateSpace && ~isempty(G_ss)
                details = [details, 'Using state-space approach for unstable plant.\n'];
                
                % Extract state-space matrices
                [A, B, C, D] = ssdata(G_ss);
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
                    
                    Wt = tf([1, w_c/1.5], [0.05, w_c*8]);  % T should be small at high frequencies
                end
                
                % Adjust weights for specific plant characteristics
                if plantInfo.isUnstable
                    % For unstable plants, increase the sensitivity weight at crossover
                    Ws = Ws * tf([1, w_c/2], [1, w_c/20]);
                    details = [details, 'Adjusted weights for unstable plant.\n'];
                    
                    % Extract unstable poles for more targeted design
                    p = eig(A);
                    unstable_idx = find(real(p) > 0);
                    
                    if ~isempty(unstable_idx)
                        % Find the most unstable pole (largest real part)
                        [max_real, idx] = max(real(p(unstable_idx)));
                        
                        % Adjust weights based on how unstable the system is
                        if max_real > 5 || length(unstable_idx) > 1
                            details = [details, 'System is highly unstable. Using more conservative design.\n'];
                            
                            % Reduce sensitivity weight to allow more control freedom
                            Ws = Ws / 2;
                            
                            % Increase control bandwidth to stabilize fast unstable modes
                            w_bc = max(w_bc, 3 * max_real);
                            Wks = tf([0.1, w_bc/2], [0.001, w_bc]);
                        end
                    end
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
                
                details = [details, sprintf('\nSensitivity weight Ws: %s\nControl effort weight Wks: %s\nRobustness weight Wt: %s\n', char(Ws), char(Wks), char(Wt))];
                
                % Setup mixed-sensitivity problem
                try
                    % Create augmented plant with weights
                    P = augw(G_ss, Ws, Wks, Wt);
                    
                    % Perform H-infinity synthesis
                    [K_hinf, CL, gamma_achieved] = hinfsyn(P, 1, 1);
                    
                    details = [details, sprintf('H-infinity synthesis successful with gamma = %.4f\n', gamma_achieved)];
                    K_temp = K_hinf;
                catch ME
                    details = [details, sprintf('H-infinity synthesis failed: %s\n', ME.message)];
                    details = [details, 'Attempting alternative approach with loop-shaping.\n'];
                    
                    % Fallback to loop-shaping approximation
                    K_temp = loopsyn(G_ss, Ws);
                    details = [details, 'Using loop-shaping approximation for H-infinity design.\n'];
                end
            else
                % Using transfer function approach
                details = [details, 'Using transfer function approach for H-infinity design.\n'];
                
                % Determine bandwidth and create weights
                try
                    w_bandwidth = getBandwidth(G, plantInfo);
                    details = [details, sprintf('Estimated plant bandwidth: %.4f rad/s\n', w_bandwidth)];
                catch
                    w_bandwidth = 1.0;
                    details = [details, 'Could not determine bandwidth. Using default value.\n'];
                end
                
                % Create mixed-sensitivity weights
                if strcmpi(robustness, 'High')
                    % Create conservative weights
                    Ws = tf([1/1.2, 0.01*w_bandwidth], [1, 0.00001*w_bandwidth]);
                    Wks = tf([1, w_bandwidth], [0.001, 10*w_bandwidth]);
                    Wt = tf([1, w_bandwidth], [0.01, 10*w_bandwidth]);
                elseif strcmpi(robustness, 'Low')
                    % Create performance-focused weights
                    Ws = tf([1/2.0, 0.05*w_bandwidth], [1, 0.0005*w_bandwidth]);
                    Wks = tf([0.1, 5*w_bandwidth], [0.01, 20*w_bandwidth]);
                    Wt = tf([1, 0.5*w_bandwidth], [0.1, 5*w_bandwidth]);
                else
                    % Create balanced weights
                    Ws = tf([1/1.5, 0.02*w_bandwidth], [1, 0.0001*w_bandwidth]);
                    Wks = tf([0.5, 2*w_bandwidth], [0.005, 15*w_bandwidth]);
                    Wt = tf([1, 0.7*w_bandwidth], [0.05, 8*w_bandwidth]);
                end
                
                % Adjust weights for unstable plants
                if plantInfo.isUnstable
                    % Make weights more aggressive for unstable plants
                    Ws = Ws * tf([1, w_bandwidth/2], [1, w_bandwidth/20]);
                    
                    % Extract unstable poles
                    unstable_poles = plantInfo.poles(real(plantInfo.poles) > 0);
                    max_real = max(real(unstable_poles));
                    
                    if max_real > 5 || length(unstable_poles) > 1
                        details = [details, 'System is highly unstable. Using more conservative design.\n'];
                        
                        % Reduce sensitivity weight
                        Ws = Ws / 2;
                        
                        % Increase control bandwidth
                        Wks = tf([0.1, 3*max_real], [0.001, 10*max_real]);
                    end
                end
                
                details = [details, sprintf('\nSensitivity weight Ws: %s\nControl effort weight Wks: %s\nRobustness weight Wt: %s\n', char(Ws), char(Wks), char(Wt))];
                
                % Try mixed-sensitivity synthesis
                try
                    % Create augmented plant
                    P = augw(G, Ws, Wks, Wt);
                    
                    % Perform H-infinity synthesis
                    [K_hinf, CL, gamma_achieved] = hinfsyn(P, 1, 1);
                    
                    details = [details, sprintf('H-infinity synthesis successful with gamma = %.4f\n', gamma_achieved)];
                    K_temp = K_hinf;
                catch ME
                    details = [details, sprintf('H-infinity synthesis failed: %s\n', ME.message)];
                    details = [details, 'Attempting alternative approach.\n'];
                    
                    % Try loop-shaping approximation
                    try
                        K_temp = loopsyn(G, Ws);
                        details = [details, 'Using loop-shaping approximation for H-infinity design.\n'];
                    catch ME2
                        details = [details, sprintf('Loop-shaping also failed: %s\n', ME2.message)];
                        details = [details, 'Using model-specific fallback approach.\n'];
                        
                        % Final fallback: create a stabilizing controller directly
                        [K_temp, ~] = createStabilizingController(G, plantInfo);
                    end
                end
            end
        else
            % For stable plants, use standard mixed-sensitivity approach
            details = [details, 'Using standard mixed-sensitivity approach for stable plant.\n'];
            
            % Determine appropriate bandwidth
            try
                w_bandwidth = getBandwidth(G, plantInfo);
                details = [details, sprintf('Estimated plant bandwidth: %.4f rad/s\n', w_bandwidth)];
            catch
                w_bandwidth = 1.0;
                details = [details, 'Could not determine bandwidth. Using default value.\n'];
            end
            
            % Create weights based on robustness setting
            if strcmpi(robustness, 'High')
                Ws = tf([1/1.2, 0.02*w_bandwidth], [1, 0.00001*w_bandwidth]);
                Wks = tf([1, w_bandwidth], [0.01, 10*w_bandwidth]);
                Wt = tf([1, w_bandwidth], [0.01, 10*w_bandwidth]);
            elseif strcmpi(robustness, 'Low')
                Ws = tf([1/2.0, 0.1*w_bandwidth], [1, 0.001*w_bandwidth]);
                Wks = tf([0.1, 2*w_bandwidth], [0.1, 20*w_bandwidth]);
                Wt = tf([1, 0.5*w_bandwidth], [0.1, 5*w_bandwidth]);
            else
                Ws = tf([1/1.5, 0.05*w_bandwidth], [1, 0.0005*w_bandwidth]);
                Wks = tf([0.5, w_bandwidth], [0.05, 15*w_bandwidth]);
                Wt = tf([1, 0.7*w_bandwidth], [0.05, 8*w_bandwidth]);
            end
            
            details = [details, sprintf('\nSensitivity weight Ws: %s\nControl effort weight Wks: %s\nRobustness weight Wt: %s\n', char(Ws), char(Wks), char(Wt))];
            
            % Try mixed-sensitivity synthesis
            try
                % Create augmented plant
                if isStateSpace && ~isempty(G_ss)
                    P = augw(G_ss, Ws, Wks, Wt);
                else
                    P = augw(G, Ws, Wks, Wt);
                end
                
                % Perform H-infinity synthesis
                [K_hinf, CL, gamma_achieved] = hinfsyn(P, 1, 1);
                
                details = [details, sprintf('H-infinity synthesis successful with gamma = %.4f\n', gamma_achieved)];
                K_temp = K_hinf;
            catch ME
                details = [details, sprintf('H-infinity synthesis failed: %s\n', ME.message)];
                details = [details, 'Using loop-shaping approximation.\n'];
                
                % Use loop-shaping as fallback
                if isStateSpace && ~isempty(G_ss)
                    K_temp = loopsyn(G_ss, Ws);
                else
                    K_temp = loopsyn(G, Ws);
                end
            end
        end
        
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
                
                % Limit to reasonable values
                Kp = min(max(abs(Kp), 0.1), 100);
                
                K = tf(Kp, 1);
                details = [details, sprintf('\nP controller with gain: Kp = %.4f', Kp)];
                
            case 'PI'
                % Extract approximate PI parameters from H-infinity controller
                try
                    % Use frequency domain approach
                    w = logspace(-4, 4, 200);
                    [mag, phase] = bode(K_temp, w);
                    mag = squeeze(mag);
                    phase = squeeze(phase);
                    
                    % Check for integral action (phase approaching -90° at low frequencies)
                    has_integrator = (phase(1) < -45);
                    
                    if has_integrator
                        % Find frequencies for integral and derivative components
                        idx_45 = find(phase > -45, 1);
                        if isempty(idx_45)
                            idx_45 = 2;
                        end
                        w_i = w(idx_45);
                        
                        % Get proportional gain at mid-frequencies
                        idx_mid = ceil(length(w)/2);
                        Kp = mag(idx_mid);
                        
                        % Calculate integral gain
                        Ki = Kp * w_i / 5;
                    else
                        % If no integral action, estimate from DC gain
                        Kp = abs(dcgain(K_temp));
                        if isnan(Kp) || isinf(Kp)
                            Kp = mag(ceil(length(w)/2));
                        end
                        
                        % Add appropriate integral action
                        Ki = Kp * w_bandwidth / 10;
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
                    
                    % Limit to reasonable values
                    Kp = min(max(abs(Kp), 0.1), 100);
                    Ki = min(max(abs(Ki), 0.01), 50);
                    
                    K = tf([Kp, Ki], [1, 0]);
                    details = [details, sprintf('\nPI controller with:\nKp = %.4f\nKi = %.4f\nTi = %.4f', Kp, Ki, Kp/Ki)];
                catch ME
                    % Fallback to simpler approach
                    Kp = abs(dcgain(K_temp));
                    if isnan(Kp) || isinf(Kp)
                        Kp = 1.0;
                    end
                    
                    % Set reasonable Ki
                    Ki = Kp * w_bandwidth / 10;
                    
                    if plantInfo.hasIntegrator
                        Ki = Ki * 0.5;
                    end
                    
                    if plantInfo.isUnstable
                        Kp = max(Kp, 1.2 * plantInfo.stabilizingGain);
                    end
                    
                    % Limit to reasonable values
                    Kp = min(max(abs(Kp), 0.1), 100);
                    Ki = min(max(abs(Ki), 0.01), 50);
                    
                    K = tf([Kp, Ki], [1, 0]);
                    details = [details, sprintf('\nSimplified PI controller with:\nKp = %.4f\nKi = %.4f', Kp, Ki)];
                end
                
            case 'PD'
                % Extract approximate PD parameters from H-infinity controller
                try
                    % Use frequency domain approach
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
                        Kp = abs(dcgain(K_temp));
                        if isnan(Kp) || isinf(Kp)
                            Kp = 1.0;
                        end
                        
                        % For plants with delay, adjust Td based on estimated delay
                        if plantInfo.hasDelay && ~isnan(plantInfo.FOPDT.L)
                            Td = 0.1 * plantInfo.FOPDT.L;
                        else
                            Td = 0.1 / w_bandwidth;  % Default
                        end
                        
                        Kd = Kp * Td;
                        details = [details, 'No clear phase lead detected. Using conservative PD parameters.\n'];
                    end
                    
                    % Adjust parameters for problematic plants
                    if plantInfo.hasRHPZeros
                        Td = Td * 0.7;  % Reduce derivative action for non-minimum phase systems
                        Kd = Kp * Td;
                        details = [details, 'Reduced derivative action for non-minimum phase plant.\n'];
                    end
                    
                    if plantInfo.isUnstable
                        Kp = max(Kp, 1.2 * plantInfo.stabilizingGain);  % Ensure stability
                        details = [details, 'Adjusted Kp to ensure stabilization.\n'];
                    end
                    
                    % Limit to reasonable values
                    Kp = min(max(abs(Kp), 0.1), 100);
                    Kd = min(max(abs(Kd), 0.01), 50);
                    
                    K = tf([Kd, Kp], [epsilon*Td, 1]);
                    details = [details, sprintf('\nPD controller with:\nKp = %.4f\nKd = %.4f\nTd = %.4f\nFilter epsilon = %.4f', Kp, Kd, Td, epsilon)];
                catch ME
                    % Fallback to simpler approach
                    Kp = abs(dcgain(K_temp));
                    if isnan(Kp) || isinf(Kp)
                        Kp = 1.0;
                    end
                    
                    Td = 0.1 / w_bandwidth;
                    Kd = Kp * Td;
                    
                    if plantInfo.hasRHPZeros
                        Kd = Kd * 0.7;
                    end
                    
                    if plantInfo.isUnstable
                        Kp = max(Kp, 1.2 * plantInfo.stabilizingGain);
                    end
                    
                    % Limit to reasonable values
                    Kp = min(max(abs(Kp), 0.1), 100);
                    Kd = min(max(abs(Kd), 0.01), 50);
                    
                    K = tf([Kd, Kp], [epsilon*Td, 1]);
                    details = [details, sprintf('\nSimplified PD controller with:\nKp = %.4f\nKd = %.4f', Kp, Kd)];
                end
                
            case 'PID'
                % Extract PID parameters from H-infinity controller frequency response
                try
                    w = logspace(-4, 4, 200);
                    [mag, phase] = bode(K_temp, w);
                    mag = squeeze(mag);
                    phase = squeeze(phase);
                    
                    % Look for integral and derivative actions
                    has_integrator = (mean(phase(1:min(5, length(phase)))) < -45);
                    has_derivative = any(phase > 10);
                    
                    if has_integrator && has_derivative
                        % Full PID behavior detected
                        details = [details, 'Full PID behavior detected in H-infinity controller.\n'];
                        
                        % Find integral and derivative frequencies
                        idx_45 = find(phase > -45, 1);
                        if isempty(idx_45)
                            idx_45 = 2;
                        end
                        w_i = w(idx_45);
                        
                        [~, idx_max] = max(phase);
                        w_d = w(idx_max);
                        
                        % Extract parameters
                        idx_mid = ceil(length(w)/2);
                        Kp = mag(idx_mid);
                        Ki = Kp * w_i / 5;
                        Kd = Kp / w_d;
                    else
                        % Partial PID behavior detected
                        details = [details, 'Incomplete PID behavior. Using approximation.\n'];
                        
                        % Extract what we can from frequency response
                        Kp = abs(dcgain(K_temp));
                        if isnan(Kp) || isinf(Kp)
                            Kp = mag(ceil(length(w)/2));
                        end
                        
                        if has_integrator
                            idx_45 = find(phase > -45, 1);
                            if isempty(idx_45)
                                idx_45 = 2;
                            end
                            w_i = w(idx_45);
                            Ki = Kp * w_i / 5;
                        else
                            Ki = Kp * w_bandwidth / 10;
                        end
                        
                        if has_derivative
                            [~, idx_max] = max(phase);
                            w_d = w(idx_max);
                            Kd = Kp / w_d;
                        else
                            Kd = Kp / w_bandwidth;
                        end
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
                        
                        % For highly unstable systems, be even more cautious with integral action
                        unstable_poles = plantInfo.poles(real(plantInfo.poles) > 0);
                        if max(real(unstable_poles)) > 5 || length(unstable_poles) > 1
                            Ki = Ki * 0.1;
                            details = [details, 'Greatly reduced integral gain for highly unstable plant.\n'];
                        end
                    end
                    
                    % Limit to reasonable values
                    Kp = min(max(abs(Kp), 0.1), 100);
                    Ki = min(max(abs(Ki), 0.01), 50);
                    Kd = min(max(abs(Kd), 0.01), 50);
                    
                    % Calculate traditional time constants
                    Ti = Kp / Ki;
                    Td = Kd / Kp;
                    
                    K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                    details = [details, sprintf('\nPID controller with:\nKp = %.4f\nKi = %.4f\nKd = %.4f\nTi = %.4f\nTd = %.4f\nFilter epsilon = %.4f', Kp, Ki, Kd, Ti, Td, epsilon)];
                catch ME
                    % Fallback to simpler PID design
                    Kp = abs(dcgain(K_temp));
                    if isnan(Kp) || isinf(Kp)
                        Kp = 1.0;
                    end
                    
                    Ki = Kp * w_bandwidth / 10;
                    Kd = Kp / w_bandwidth;
                    
                    if plantInfo.hasIntegrator
                        Ki = Ki * 0.5;
                    end
                    
                    if plantInfo.hasRHPZeros
                        Kd = Kd * 0.7;
                    end
                    
                    if plantInfo.isUnstable
                        Kp = max(Kp, 1.2 * plantInfo.stabilizingGain);
                        
                        unstable_poles = plantInfo.poles(real(plantInfo.poles) > 0);
                        if max(real(unstable_poles)) > 5 || length(unstable_poles) > 1
                            Ki = Ki * 0.1;
                        end
                    end
                    
                    % Limit to reasonable values
                    Kp = min(max(abs(Kp), 0.1), 100);
                    Ki = min(max(abs(Ki), 0.01), 50);
                    Kd = min(max(abs(Kd), 0.01), 50);
                    
                    K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                    details = [details, sprintf('\nSimplified PID controller with:\nKp = %.4f\nKi = %.4f\nKd = %.4f', Kp, Ki, Kd)];
                end
                
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
                    details = [details, sprintf('\nClosed-loop performance:\nGain Margin: %.2f dB\nPhase Margin: %.2f degrees', 20*log10(Gm), Pm)];
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
            % For unstable plants, use a stabilizing controller
            [K, fallback_details] = createStabilizingController(G, plantInfo);
            details = [details, fallback_details];
        else
            % For stable plants, use a conservative controller
            switch structure
                case 'P'
                    Kp = 0.5;
                    if ~isnan(plantInfo.dcGain) && ~isinf(plantInfo.dcGain) && plantInfo.dcGain ~= 0
                        Kp = 0.5 / abs(plantInfo.dcGain);
                    end
                    K = tf(Kp, 1);
                    details = [details, sprintf('Conservative P controller with Kp = %.4f', Kp)];
                    
                case 'PI'
                    Kp = 0.5;
                    if ~isnan(plantInfo.dcGain) && ~isinf(plantInfo.dcGain) && plantInfo.dcGain ~= 0
                        Kp = 0.5 / abs(plantInfo.dcGain);
                    end
                    Ki = Kp * 0.1;
                    K = tf([Kp, Ki], [1, 0]);
                    details = [details, sprintf('Conservative PI controller with Kp = %.4f, Ki = %.4f', Kp, Ki)];
                    
                case 'PD'
                    Kp = 0.5;
                    if ~isnan(plantInfo.dcGain) && ~isinf(plantInfo.dcGain) && plantInfo.dcGain ~= 0
                        Kp = 0.5 / abs(plantInfo.dcGain);
                    end
                    Kd = Kp * 0.1;
                    K = tf([Kd, Kp], [epsilon*Kd, 1]);
                    details = [details, sprintf('Conservative PD controller with Kp = %.4f, Kd = %.4f', Kp, Kd)];
                    
                case 'PID'
                    Kp = 0.5;
                    if ~isnan(plantInfo.dcGain) && ~isinf(plantInfo.dcGain) && plantInfo.dcGain ~= 0
                        Kp = 0.5 / abs(plantInfo.dcGain);
                    end
                    Ki = Kp * 0.1;
                    Kd = Kp * 0.1;
                    K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                    details = [details, sprintf('Conservative PID controller with Kp = %.4f, Ki = %.4f, Kd = %.4f', Kp, Ki, Kd)];
                    
                otherwise
                    K = tf(0.5, 1);
                    details = [details, 'Using conservative P controller as fallback.'];
            end
        end
    end
end

% Helper function to create a stabilizing controller for unstable plants
function [K, details] = createStabilizingController(G, plantInfo)
    % Create a stabilizing controller based on plant info
    details = '';
    
    % Get unstable poles
    p = plantInfo.poles;
    unstable_poles = p(real(p) > 0);
    
    if isempty(unstable_poles)
        % If plant is actually stable, use a simple controller
        K = tf(0.5, 1);
        details = 'Using simple controller for stable plant.';
        return;
    end
    
    % For unstable plants, create a pole-zero cancellation controller
    K_num = 1;
    K_den = 1;
    
    for i = 1:length(unstable_poles)
        pole_i = unstable_poles(i);
        
        if imag(pole_i) ~= 0
            % Skip complex conjugate pairs, we'll handle them together
            if i < length(unstable_poles) && abs(pole_i - conj(unstable_poles(i+1))) < 1e-6
                continue;
            end
            
            % For complex poles, create zeros that cancel the pole pair
            if imag(pole_i) > 0
                real_part = real(pole_i);
                imag_part = imag(pole_i);
                
                % Zero to cancel pole pair: s^2 - 2*re(p)*s + |p|^2
                num_term = [1, -2*real_part, real_part^2 + imag_part^2];
                
                % Stable poles with the same natural frequency but more damping
                den_term = [1, 2*abs(real_part), real_part^2 + imag_part^2];
                
                K_num = conv(K_num, num_term);
                K_den = conv(K_den, den_term);
            end
        else
            % Real unstable pole
            K_num = conv(K_num, [1, -pole_i]);
            K_den = conv(K_den, [1, abs(pole_i)]);
        end
    end
    
    % Calculate appropriate gain
    max_real_part = max(real(unstable_poles));
    K_gain = max(1.0, max_real_part);
    
    K = tf(K_gain * K_num, K_den);
    details = sprintf('Created stabilizing controller for unstable plant with gain = %.4f', K_gain);
    
    % Test if controller stabilizes the plant
    try
        CL = feedback(G*K, 1);
        cl_poles = pole(CL);
        
        if any(real(cl_poles) > 0)
            % Try reducing gain
            for scale = [0.5, 0.2, 0.1, 0.05, 0.01]
                K_test = tf(K_gain * scale * K_num, K_den);
                CL_test = feedback(G * K_test, 1);
                
                if all(real(pole(CL_test)) < 0)
                    K = K_test;
                    details = [details, sprintf('\nGain scaled by %.2f to achieve stability.', scale)];
                    break;
                end
            end
        end
    catch
        % If feedback analysis fails, keep original controller
    end
end

% Function to get a formatted string with plant information
function infoStr = getPlantInfoString(plantInfo)
    % Initialize output string
    infoStr = '';
    
    % Add stability information
    if plantInfo.isUnstable
        infoStr = [infoStr, 'Unstable, '];
    else
        infoStr = [infoStr, 'Stable, '];
    end
    
    % Add integrator information
    if plantInfo.hasIntegrator
        infoStr = [infoStr, 'Has integrator, '];
    end
    
    % Add RHP zeros information
    if plantInfo.hasRHPZeros
        infoStr = [infoStr, 'Non-minimum phase, '];
    end
    
    % Add delay information
    if plantInfo.hasDelay
        infoStr = [infoStr, 'Has delay, '];
    end
    
    % Add order information
    if plantInfo.isHighOrder
        infoStr = [infoStr, 'High-order, '];
    else
        infoStr = [infoStr, 'Low-order, '];
    end
    
    % Add DC gain if available
    if ~isnan(plantInfo.dcGain) && ~isinf(plantInfo.dcGain)
        infoStr = [infoStr, sprintf('DC gain=%.3g', plantInfo.dcGain)];
    else
        infoStr = [infoStr, 'Infinite DC gain'];
    end
end