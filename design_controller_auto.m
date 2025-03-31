function [K, details] = design_controller_auto(G, method, structure, options)
% DESIGN_CONTROLLER_AUTO Automatic controller design with various methods
% 
% Inputs:
%   G         - Plant transfer function
%   method    - Design method (string)
%   structure - Controller structure: 'P', 'PI', 'PD', 'PID'
%   options   - Structure with optional parameters:
%     .epsilon     - Filter parameter for D-term (default: 0.1)
%     .phaseMargin - Desired phase margin in degrees (default: 45)
%     .bandwidth   - Desired bandwidth in rad/s (default: 1)
%     .settlingTime - Desired settling time in s (default: 5)
%     .robustness  - Robustness level: 'Low', 'Medium', 'High' (default: 'Medium')
%     .overshoot   - Desired overshoot in % (default: 10)
%     .goal        - Optimization goal: 'Tracking', 'Disturbance Rejection', 'Robustness' (default: 'Tracking')
%
% Outputs:
%   K        - Designed controller as transfer function
%   details  - Text description with design details

    % Default values for missing options
    if nargin < 4
        options = struct();
    end
    
    if ~isfield(options, 'epsilon')
        options.epsilon = 0.1;
    end
    
    if ~isfield(options, 'phaseMargin')
        options.phaseMargin = 45;
    end
    
    if ~isfield(options, 'bandwidth')
        options.bandwidth = 1;
    end
    
    if ~isfield(options, 'settlingTime')
        options.settlingTime = 5;
    end
    
    if ~isfield(options, 'robustness')
        options.robustness = 'Medium';
    end
    
    if ~isfield(options, 'overshoot')
        options.overshoot = 10;
    end
    
    if ~isfield(options, 'goal')
        options.goal = 'Tracking';
    end
    
    % NEW: Pre-analyze the plant to determine characteristics
    plantInfo = analyzePlant(G);
    
    % Method selection based on plant analysis
    if plantInfo.isUnstable && ~contains(method, {'H-infinity', 'Loop-Shaping'})
        % For unstable plants, display warning if not using robust methods
        warning(['Plant is unstable. Consider using H-infinity or Loop-Shaping methods. ', ...
                'Attempting to modify %s method to handle instability.'], method);
        
        % For unstable systems, use more robust tuning parameters
        options.robustness = 'High';
    end
    
    % Select design method
    try
        switch method
            case 'Ziegler-Nichols (Oscillation)'
                [K, details] = designZieglerNicholsOscillation(G, structure, options.epsilon, plantInfo);
            case 'Ziegler-Nichols (Step)'
                [K, details] = designZieglerNicholsStep(G, structure, options.epsilon, plantInfo);
            case 'Aström'
                [K, details] = designAstrom(G, structure, options.epsilon, plantInfo);
            case 'CHR'
                [K, details] = designCHR(G, structure, options.epsilon, plantInfo);
            case 'Cohen-Coon'
                [K, details] = designCohenCoon(G, structure, options.epsilon, plantInfo);
            case 'Loop-Shaping'
                [K, details] = designLoopShaping(G, structure, options.phaseMargin, options.bandwidth, options.epsilon, plantInfo);
            case 'IMC'
                [K, details] = designIMC(G, structure, options.settlingTime, options.epsilon, plantInfo);
            case 'MIGO'
                [K, details] = designMIGO(G, structure, options.robustness, options.epsilon, plantInfo);
            case 'H-infinity'
                [K, details] = designHInfinity(G, structure, options.robustness, options.epsilon, plantInfo);
            case 'LQG (Linear-Quadratic-Gaussian)'
                [K, details] = designLQG(G, structure, options.bandwidth, options.robustness, options.epsilon, plantInfo);
            otherwise
                error('Unknown design method: %s', method);
        end
    catch ME
        % Enhanced error handling with fallback controller design
        warning('Error in %s design: %s\nUsing robust fallback design.', method, ME.message);
        [K, details] = designFallbackController(G, structure, options.epsilon, plantInfo);
        details = sprintf('Original design failed: %s\n\n%s', ME.message, details);
    end
    
    % Evaluate controller quality if requested
    try
        score = evaluateController(K, G, options.goal, options.phaseMargin, options.overshoot, options.settlingTime, options.bandwidth);
        details = [details, sprintf('\n\nController Score: %.2f/100', score)];
    catch ME
        disp(['Warning: Could not evaluate controller quality: ', ME.message]);
    end
end

function plantInfo = analyzePlant(G)
    % ANALYZEPLANT Analyze plant characteristics to guide controller design
    % Outputs a structure with plant information
    
    plantInfo = struct();
    
    % Get poles and zeros
    plantInfo.poles = pole(G);
    try
        plantInfo.zeros = zero(G);
    catch
        plantInfo.zeros = [];
    end
    
    % Check stability
    plantInfo.isUnstable = any(real(plantInfo.poles) > 0);
    
    % Check for integrators (poles at the origin)
    plantInfo.hasIntegrator = any(abs(plantInfo.poles) < 1e-6);
    
    % Check for non-minimum phase zeros (RHP zeros)
    plantInfo.hasRHPZeros = any(real(plantInfo.zeros) > 0);
    
    % Check for time delay
    [num, den] = tfdata(G, 'v');
    % Pade approximations typically have alternating sign coefficients
    if length(num) > 1 && all(sign(num(1:2:end)) ~= sign(num(2:2:end)))
        plantInfo.hasDelay = true;
    else
        plantInfo.hasDelay = false;
    end
    
    % Determine plant DC gain
    try
        plantInfo.dcGain = dcgain(G);
    catch
        % For plants with pure integrators
        plantInfo.dcGain = Inf;
    end
    
    % Check if high order (more than 2 states)
    plantInfo.isHighOrder = length(plantInfo.poles) > 2;
    
    % Get step response characteristics if plant is stable
    if ~plantInfo.isUnstable
        try
            t = linspace(0, 100, 1000);
            [y, t] = step(G, t);
            info = stepinfo(y, t);
            plantInfo.stepInfo = info;
            plantInfo.stepResponse = struct('time', t, 'response', y);
            
            % Compute approximate first-order + delay model
            [L, T, K] = estimateFirstOrderPlusDelay(t, y);
            plantInfo.FOPDT = struct('L', L, 'T', T, 'K', K);
        catch
            plantInfo.stepInfo = [];
            plantInfo.stepResponse = [];
            plantInfo.FOPDT = struct('L', NaN, 'T', NaN, 'K', NaN);
        end
    else
        plantInfo.stepInfo = [];
        plantInfo.stepResponse = [];
        plantInfo.FOPDT = struct('L', NaN, 'T', NaN, 'K', NaN);
    end
    
    % For unstable plants, estimate stabilizing feedback
    if plantInfo.isUnstable
        plantInfo.stabilizingGain = estimateStabilizingGain(G);
    else
        plantInfo.stabilizingGain = 0;
    end
end

function K = preStabilize(G, plantInfo)
    % PRESTABILIZE Create a simple stabilizing controller for unstable plants
    
    if ~plantInfo.isUnstable
        K = tf(1, 1); % Unity controller for already stable plants
        return;
    end
    
    % Start with a proportional controller
    K = tf(-plantInfo.stabilizingGain, 1);
    
    % Add small derivative action if high-frequency dynamics present
    if plantInfo.isHighOrder
        Td = 0.05;
        epsilon = 0.1;
        K = tf([Td, 1], [epsilon*Td, 1]) * K;
    end
    
    % Check if our stabilization worked
    try
        closedloop = feedback(G*K, 1);
        poles = pole(closedloop);
        if any(real(poles) > 0)
            % If still unstable, try a different approach with PD control
            Kp = -plantInfo.stabilizingGain * 1.5;
            Td = 0.1;
            epsilon = 0.1;
            K = tf([Td*Kp, Kp], [epsilon*Td, 1]);
        end
    catch
        % If feedback fails, use a more conservative stabilizing controller
        K = tf(-2*plantInfo.stabilizingGain, [0.1, 1]);
    end
end

function stabilizingGain = estimateStabilizingGain(G)
    % ESTIMATESTABILIZINGGAIN Estimate a stabilizing feedback gain for unstable plants
    
    p = pole(G);
    
    % Find the most unstable pole
    [maxRealPart, idx] = max(real(p));
    
    if maxRealPart <= 0
        % Plant is stable
        stabilizingGain = 0;
        return;
    end
    
    % For a dominant real unstable pole, we need a proportional gain
    % that moves it to the LHP. Use a simple mapping formula.
    if imag(p(idx)) == 0
        % Pure real pole - simple proportional stabilization
        stabilizingGain = maxRealPart * 1.5;
    else
        % Complex pole - use a more conservative gain
        stabilizingGain = maxRealPart * 2;
    end
end

function [L, T, K] = estimateFirstOrderPlusDelay(t, y)
    % ESTIMATEFIRSTORDERPLUSDELAY Estimate a FOPDT model from step response
    
    % Get steady-state gain
    y_final = y(end);
    K = y_final;
    
    if abs(K) < 1e-6
        % No steady-state response
        L = 0.1;
        T = 1.0;
        return;
    end
    
    % Normalize the response
    y_norm = y / y_final;
    
    % Find 63.2% response point for time constant
    idx_63 = find(y_norm >= 0.632, 1);
    
    if isempty(idx_63)
        % Use 95% of response time as approximation
        idx_95 = find(y_norm >= 0.95, 1);
        if isempty(idx_95)
            T = t(end) / 3;
        else
            T = t(idx_95) / 3;
        end
    else
        T = t(idx_63);
    end
    
    % Find 10% response for delay estimation
    idx_10 = find(y_norm >= 0.1, 1);
    
    if isempty(idx_10)
        L = T * 0.1; % Default approximation
    else
        L = t(idx_10);
    end
    
    % Ensure reasonable values
    L = max(0.01, min(L, T)); % Bound delay to be less than time constant
    T = max(0.1, T);
end

function [K, details] = designZieglerNicholsOscillation(G, structure, epsilon, plantInfo)
    % Improved Ziegler-Nichols oscillation method with pre-stabilization for unstable plants
    
    % Add plant information to details
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % For unstable plants, pre-stabilize before Z-N tuning
    if plantInfo.isUnstable
        K_stab = preStabilize(G, plantInfo);
        details = [details, sprintf('Plant is unstable. Using pre-stabilization.\n')];
        G_stab = feedback(G, K_stab);
    else
        G_stab = G;
    end
    
    % Determine critical gain and period
    k_krit = 0;
    T_krit = 0;
    
    try
        % Use frequency domain approach for critical point
        w = logspace(-2, 4, 1000);
        [mag, phase] = bode(G_stab, w);
        mag = squeeze(mag);
        phase = squeeze(phase);
        
        % Find where phase crosses -180 degrees
        phase_crossings = [];
        for i = 1:length(phase)-1
            if (phase(i) > -180 && phase(i+1) < -180) || (phase(i) < -180 && phase(i+1) > -180)
                phase_crossings(end+1) = i;
            end
        end
        
        if ~isempty(phase_crossings)
            % Linear interpolation to find precise crossing frequency
            for crossing = phase_crossings
                w1 = w(crossing);
                w2 = w(crossing+1);
                phase1 = phase(crossing);
                phase2 = phase(crossing+1);
                
                w_cross = w1 + (w2 - w1) * (-180 - phase1) / (phase2 - phase1);
                
                % Get magnitude at this frequency using interpolation
                mag1 = mag(crossing);
                mag2 = mag(crossing+1);
                mag_cross = mag1 + (mag2 - mag1) * (w_cross - w1) / (w2 - w1);
                
                % Critical gain is reciprocal of magnitude at -180 degrees
                k_krit_candidate = 1 / mag_cross;
                T_krit_candidate = 2*pi / w_cross;
                
                % Use the first valid crossing point
                if k_krit == 0 || k_krit_candidate < k_krit
                    k_krit = k_krit_candidate;
                    T_krit = T_krit_candidate;
                end
            end
        end
        
        % If frequency domain approach failed, use time-domain approach
        if k_krit == 0
            % Start with a small gain and increase it iteratively
            found = false;
            step_size = 0.1;
            max_iterations = 2000;
            
            for iteration = 1:max_iterations
                k = iteration * step_size;
                
                % Test stability with current k
                closed_loop = feedback(G_stab*k, 1);
                poles = pole(closed_loop);
                
                % Check for poles on the imaginary axis
                realParts = real(poles);
                imagParts = imag(poles);
                
                close_to_imag_axis = find(abs(realParts) < 0.001 & imagParts ~= 0);
                
                if ~isempty(close_to_imag_axis)
                    % Found poles near stability boundary
                    k_krit = k;
                    
                    % Calculate period
                    idx = find(abs(realParts) < 0.001 & imagParts > 0, 1);
                    if ~isempty(idx)
                        omega = imagParts(idx);
                        T_krit = 2*pi/omega;
                        found = true;
                        break;
                    end
                elseif any(realParts > 0)
                    % System became unstable
                    
                    % Reduce k incrementally to find boundary
                    for back_step = 1:10
                        k_test = k - back_step * (step_size/10);
                        
                        closed_loop = feedback(G_stab*k_test, 1);
                        poles = pole(closed_loop);
                        realParts = real(poles);
                        imagParts = imag(poles);
                        
                        if all(realParts < 0) && any(abs(realParts) < 0.01 & imagParts ~= 0)
                            k_krit = k_test;
                            
                            idx = find(abs(realParts) < 0.01 & imagParts > 0, 1);
                            if ~isempty(idx)
                                omega = imagParts(idx);
                                T_krit = 2*pi/omega;
                                found = true;
                                break;
                            end
                        end
                    end
                    
                    if found
                        break;
                    else
                        % Assume we've passed the critical point
                        k_krit = k - step_size;
                        
                        closed_loop = feedback(G_stab*k_krit, 1);
                        poles = pole(closed_loop);
                        imagParts = imag(poles(abs(real(poles)) < 0.1));
                        
                        if ~isempty(imagParts) && any(imagParts > 0)
                            omega = max(imagParts(imagParts > 0));
                            T_krit = 2*pi/omega;
                            found = true;
                            break;
                        end
                    end
                    
                    break;
                end
            end
        end
        
        if k_krit == 0
            % Fallback estimation based on plant characteristics
            if ~isnan(plantInfo.FOPDT.L) && ~isnan(plantInfo.FOPDT.T)
                % Approximate critical values based on FOPDT model
                k_krit = (plantInfo.FOPDT.T / plantInfo.FOPDT.K / plantInfo.FOPDT.L) * pi / 2;
                T_krit = 4 * plantInfo.FOPDT.L;
                details = [details, sprintf('Using estimated critical parameters from FOPDT model.\n')];
            else
                error('Could not determine critical gain. The plant may not be suitable for the Ziegler-Nichols oscillation method.');
            end
        end
        
    catch ME
        % If all else fails, use approximation based on poles
        if plantInfo.isHighOrder
            % For higher-order systems, use conservative approximation
            k_krit = 2.0;
            T_krit = 1.0;
            details = [details, sprintf('Critical parameters estimated due to error: %s\n', ME.message)];
        else
            error('Error applying Ziegler-Nichols oscillation method: %s', ME.message);
        end
    end
    
    % Controller parameters according to Ziegler-Nichols table
    details = [details, sprintf('k_krit = %.4f\nT_krit = %.4f s\n', k_krit, T_krit)];
    
    switch structure
        case 'P'
            Kp = 0.5 * k_krit;
            K = tf(Kp, 1);
        case 'PI'
            Kp = 0.45 * k_krit;
            Ti = 0.85 * T_krit;
            K = tf([Kp, Kp/Ti], [1, 0]);
        case 'PD'
            Kp = 0.5 * k_krit;
            Td = 0.12 * T_krit;
            K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
        case 'PID'
            % Use modified parameters for better robustness
            Kp = 0.6 * k_krit;
            Ti = 0.5 * T_krit;
            Td = 0.125 * T_krit;
            
            % For plants with RHP zeros, reduce derivative action
            if plantInfo.hasRHPZeros
                Td = Td * 0.5;
                details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
            end
            
            K = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
        otherwise
            error('Unsupported controller structure for Ziegler-Nichols oscillation method');
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
    
    % Verify controller stability
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

function [K, details] = designCHR(G, structure, epsilon, plantInfo)
    % Enhanced CHR (Chien-Hrones-Reswick) method for 0% overshoot with improved stability handling
    
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
            error('System does not have a finite DC gain. Not suitable for CHR method.');
        end
    end
    
    % Normalize the response
    y_norm = y / y_final;
    
    % Use FOPDT parameters if available, otherwise estimate them
    if ~isnan(plantInfo.FOPDT.L) && ~isnan(plantInfo.FOPDT.T)
        L = plantInfo.FOPDT.L;
        T = plantInfo.FOPDT.T;
        Ks = plantInfo.FOPDT.K;
    else
        try
            % Determine the tangent at the inflection point
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
            % Fallback method if tangent method fails
            idx_63 = find(y_norm >= 0.632, 1);
            if isempty(idx_63)
                error('Could not determine time constant. The plant may not be suitable for the CHR method.');
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
        end
    end
    
    % Calculate parameter ratio for tuning adjustment
    ratio = L/T;
    
    details = [details, sprintf('Ks = %.4f\nL = %.4f s\nT = %.4f s\nL/T ratio = %.4f\n', Ks, L, T, ratio)];
    
    % Apply adjustment factor based on plant characteristics
    adjustmentFactor = 1.0;
    
    if ratio > 0.5
        % High delay-to-time-constant ratio - use more conservative settings
        adjustmentFactor = 0.8;
        details = [details, 'High L/T ratio: Using more conservative tuning.\n'];
    elseif plantInfo.hasRHPZeros
        % Non-minimum phase - use more conservative settings
        adjustmentFactor = 0.7;
        details = [details, 'Non-minimum phase behavior: Using more conservative tuning.\n'];
    elseif plantInfo.hasIntegrator
        % Plant with integrator - reduce integral action
        adjustmentFactor = 0.9;
        details = [details, 'Plant with integrator: Adjusting parameters.\n'];
    end
    
    % Parameters according to CHR table for setpoint regulation (0% overshoot)
    switch structure
        case 'P'
            Kp = adjustmentFactor * 0.3 * T / (Ks * L);
            K = tf(Kp, 1);
        case 'PI'
            Kp = adjustmentFactor * 0.35 * T / (Ks * L);
            Ti = 1.2 * T;
            K = tf([Kp, Kp/Ti], [1, 0]);
        case 'PD'
            Kp = adjustmentFactor * 0.5 * T / (Ks * L);
            Td = 0.3 * L;
            
            % For plants with RHP zeros, reduce derivative action
            if plantInfo.hasRHPZeros
                Td = Td * 0.5;
                details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
            end
            
            K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
        case 'PID'
            Kp = adjustmentFactor * 0.6 * T / (Ks * L);
            Ti = T;
            Td = 0.5 * L;
            
            % For plants with RHP zeros, reduce derivative action
            if plantInfo.hasRHPZeros
                Td = Td * 0.5;
                details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
            end
            
            K = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
        otherwise
            error('Unsupported controller structure for CHR method');
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

function [K, details] = designCohenCoon(G, structure, epsilon, plantInfo)
    % Enhanced Cohen-Coon method with improved stability handling and parameter selection
    
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
            error('System does not have a finite DC gain. Not suitable for Cohen-Coon method.');
        end
    end
    
    % Normalize the response
    y_norm = y / y_final;
    
    % Use FOPDT parameters if available, otherwise estimate them
    if ~isnan(plantInfo.FOPDT.L) && ~isnan(plantInfo.FOPDT.T)
        L = plantInfo.FOPDT.L;
        T = plantInfo.FOPDT.T;
        Ks = plantInfo.FOPDT.K;
    else
        try
            % Determine the tangent at the inflection point
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
            % Fallback method if tangent method fails
            idx_63 = find(y_norm >= 0.632, 1);
            if isempty(idx_63)
                error('Could not determine time constant. The plant may not be suitable for the Cohen-Coon method.');
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
        end
    end
    
    % Parameter tau = T/L, helpful for Cohen-Coon formulas
    tau = T/L;
    
    details = [details, sprintf('Ks = %.4f\nL = %.4f s\nT = %.4f s\ntau = %.4f\n', Ks, L, T, tau)];
    
    % Apply adjustment factor based on plant characteristics
    adjustmentFactor = 1.0;
    
    if tau < 3
        % Low tau - reduce controller gain for stability
        adjustmentFactor = 0.8;
        details = [details, 'Low tau value: Reducing controller gain for stability.\n'];
    elseif plantInfo.hasRHPZeros
        % Non-minimum phase - use more conservative settings
        adjustmentFactor = 0.7;
        details = [details, 'Non-minimum phase behavior: Using more conservative tuning.\n'];
    elseif plantInfo.hasDelay
        adjustmentFactor = 0.85;
        details = [details, 'Detected significant delay: Adjusting parameters.\n'];
    elseif plantInfo.isHighOrder
        adjustmentFactor = 0.9;
        details = [details, 'High-order system: Adjusting for robustness.\n'];
    end
    
    % Improved Cohen-Coon formulas with adjustments
    switch structure
        case 'P'
            Kp = adjustmentFactor * (1/Ks) * (1 + (1/3)*tau);
            K = tf(Kp, 1);
        case 'PI'
            Kp = adjustmentFactor * (1/Ks) * (0.9 + (1/12)*tau);
            Ti = L * ((30 + 3*tau)/(9 + 20*tau));
            
            % For plants with integrator, adjust Ti
            if plantInfo.hasIntegrator
                Ti = Ti * 1.5;
                details = [details, 'Adjusted Ti for plant with integrator.\n'];
            end
            
            K = tf([Kp, Kp/Ti], [1, 0]);
        case 'PD'
            Kp = adjustmentFactor * (1/Ks) * (1.25 * (1 + (1/6)*tau));
            Td = L * ((6 - 2*tau)/(22 + 3*tau));
            
            % For plants with RHP zeros, reduce derivative action
            if plantInfo.hasRHPZeros
                Td = Td * 0.5;
                details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
            end
            
            % Ensure Td is positive
            Td = max(Td, 0.05*L);
            
            K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
        case 'PID'
            Kp = adjustmentFactor * (1/Ks) * (1.35 + (1/4)*tau);
            Ti = L * ((32 + 6*tau)/(13 + 8*tau));
            Td = L * (4/(11 + 2*tau));
            
            % For plants with RHP zeros, reduce derivative action
            if plantInfo.hasRHPZeros
                Td = Td * 0.5;
                details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
            end
            
            % For plants with integrator, adjust Ti
            if plantInfo.hasIntegrator
                Ti = Ti * 1.5;
                details = [details, 'Adjusted Ti for plant with integrator.\n'];
            end
            
            K = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
        otherwise
            error('Unsupported controller structure for Cohen-Coon method');
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

function [K, details] = designMIGO(G, structure, robustness, epsilon, plantInfo)
    % Enhanced MIGO (M-constrained Integral Gain Optimization) method
    
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % Set M-constraint based on robustness level
    switch robustness
        case 'Low'
            M_s = 2.0;  % Lower robustness, higher performance
        case 'Medium'
            M_s = 1.5;  % Medium robustness
        case 'High'
            M_s = 1.2;  % Higher robustness, lower performance
        otherwise
            M_s = 1.5;  % Default value
    end
    
    details = [details, sprintf('Robustness level: %s\nM_s constraint: %.2f\n', robustness, M_s)];
    
    % For unstable plants, pre-stabilize first
    if plantInfo.isUnstable
        details = [details, 'Plant is unstable. Using pre-stabilization approach.\n'];
        
        K_stab = preStabilize(G, plantInfo);
        G_stab = feedback(G, K_stab);
        
        % Continue MIGO design with stabilized plant
        G_for_design = G_stab;
    else
        G_for_design = G;
    end
    
    % Implementation based on controller structure
    switch structure
        case 'P'
            % For P-controller: find maximum Kp under M_s constraint
            Kp = findOptimalKp(G_for_design, M_s);
            K = tf(Kp, 1);
            details = [details, sprintf('Optimal P-controller gain: Kp = %.4f\n', Kp)];
            
        case 'PI'
            % For PI-controller: optimize Kp and Ki under M_s constraint
            [Kp, Ki] = findOptimalPI(G_for_design, M_s, plantInfo);
            K = tf([Kp, Ki], [1, 0]);
            Ti = Kp/Ki;
            details = [details, sprintf('Optimal PI parameters:\nKp = %.4f\nKi = %.4f\nTi = %.4f\n', Kp, Ki, Ti)];
            
        case 'PD'
            % For PD-controller: optimize Kp and Kd
            [Kp, Kd] = findOptimalPD(G_for_design, M_s, epsilon, plantInfo);
            K = tf([Kd, Kp], [epsilon*Kd, 1]);
            Td = Kd/Kp;
            details = [details, sprintf('Optimal PD parameters:\nKp = %.4f\nKd = %.4f\nTd = %.4f\n', Kp, Kd, Td)];
            
        case 'PID'
            % For PID-controller: optimize Kp, Ki, and Kd under M_s constraint
            [Kp, Ki, Kd] = findOptimalPID(G_for_design, M_s, epsilon, plantInfo);
            K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
            Ti = Kp/Ki;
            Td = Kd/Kp;
            details = [details, sprintf('Optimal PID parameters:\nKp = %.4f\nKi = %.4f\nKd = %.4f\nTi = %.4f\nTd = %.4f\n', ...
                      Kp, Ki, Kd, Ti, Td)];
            
        otherwise
            error('MIGO is only implemented for P, PI, PD, and PID controllers');
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
        end
    catch
        details = [details, 'Could not verify closed-loop stability.\n'];
    end
end

function Kp = findOptimalKp(G, M_s)
    % Find the maximum P-controller gain that satisfies the M_s constraint
    
    % Define search range
    Kp_min = 0.01;
    Kp_max = 100;
    Kp_values = logspace(log10(Kp_min), log10(Kp_max), 50);
    
    Kp_best = Kp_min;
    
    for i = 1:length(Kp_values)
        Kp = Kp_values(i);
        K_test = tf(Kp, 1);
        L = G * K_test;
        S = feedback(1, L);
        
        try
            % Calculate maximum sensitivity
            [peakgain, ~, ~] = getPeakGain(S);
            
            if peakgain <= M_s && Kp > Kp_best
                Kp_best = Kp;
            end
        catch
            % Skip if peak gain calculation fails
            continue;
        end
    end
    
    Kp = Kp_best;
end

function [Kp, Ki] = findOptimalPI(G, M_s, plantInfo)
    % Find optimal PI controller parameters under M_s constraint
    
    % Grid search approach with refinement
    % Start with coarse grid
    Kp_values = logspace(-2, 2, 15);
    Ti_values = logspace(-2, 2, 15);
    
    best_Kp = 0.1;
    best_Ti = 1.0;
    best_score = -Inf;
    
    % First pass: coarse grid search
    for i = 1:length(Kp_values)
        for j = 1:length(Ti_values)
            Kp = Kp_values(i);
            Ti = Ti_values(j);
            Ki = Kp / Ti;
            
            K_test = tf([Kp, Ki], [1, 0]);
            L = G * K_test;
            S = feedback(1, L);
            
            try
                % Calculate maximum sensitivity
                [peakgain, ~, ~] = getPeakGain(S);
                
                if peakgain <= M_s
                    % Score based on combined performance metrics
                    [~, Pm] = margin(L);
                    
                    % Calculate step response metrics if possible
                    try
                        T = feedback(G * K_test, 1);
                        step_info = stepinfo(T);
                        settling_time = step_info.SettlingTime;
                        overshoot = step_info.Overshoot;
                        
                        % Calculate score based on multiple criteria
                        % - Higher Ki for better tracking
                        % - Higher phase margin for stability
                        % - Lower overshoot and settling time for better performance
                        score = Ki * (Pm/100) * (100/(1 + overshoot)) * (10/(1 + settling_time));
                    catch
                        % If step response fails, use simpler scoring
                        score = Ki * (Pm/100);
                    end
                    
                    if score > best_score
                        best_score = score;
                        best_Kp = Kp;
                        best_Ti = Ti;
                    end
                end
            catch
                % Skip if analysis fails
                continue;
            end
        end
    end
    
    % Second pass: refined search around best point
    Kp_min = best_Kp * 0.5;
    Kp_max = best_Kp * 2.0;
    Ti_min = best_Ti * 0.5;
    Ti_max = best_Ti * 2.0;
    
    Kp_refined = linspace(Kp_min, Kp_max, 10);
    Ti_refined = linspace(Ti_min, Ti_max, 10);
    
    for i = 1:length(Kp_refined)
        for j = 1:length(Ti_refined)
            Kp = Kp_refined(i);
            Ti = Ti_refined(j);
            Ki = Kp / Ti;
            
            K_test = tf([Kp, Ki], [1, 0]);
            L = G * K_test;
            S = feedback(1, L);
            
            try
                % Calculate maximum sensitivity
                [peakgain, ~, ~] = getPeakGain(S);
                
                if peakgain <= M_s
                    % Score based on phase margin and integral gain
                    [~, Pm] = margin(L);
                    score = Ki * (Pm/100);
                    
                    if score > best_score
                        best_score = score;
                        best_Kp = Kp;
                        best_Ti = Ti;
                    end
                end
            catch
                % Skip if analysis fails
                continue;
            end
        end
    end
    
    % Adapt to plant characteristics
    if plantInfo.hasIntegrator
        % For plants with integrator, reduce integral action
        best_Ti = best_Ti * 2;
    end
    
    if plantInfo.hasRHPZeros
        % For non-minimum phase plants, reduce controller gain
        best_Kp = best_Kp * 0.8;
    end
    
    Kp = best_Kp;
    Ki = Kp / best_Ti;
end
    
function [Kp, Kd] = findOptimalPD(G, M_s, epsilon, plantInfo)
    % Find optimal PD controller parameters under M_s constraint
    
    % Grid search approach
    Kp_values = logspace(-2, 2, 15);
    Td_values = logspace(-2, 1, 15);
    
    best_Kp = 0.1;
    best_Td = 0.1;
    best_Pm = 0;
    
    for i = 1:length(Kp_values)
        for j = 1:length(Td_values)
            Kp = Kp_values(i);
            Td = Td_values(j);
            Kd = Kp * Td;
            
            % Fix: Use epsilon*Td instead of epsilon*Kd in the denominator
            K_test = tf([Kd, Kp], [epsilon*Td, 1]);
            L = G * K_test;
            
            try
                % Calculate maximum sensitivity and phase margin
                S = feedback(1, L);
                [peakgain, ~, ~] = getPeakGain(S);
                [~, Pm] = margin(L);
                
                if peakgain <= M_s && Pm > best_Pm
                    best_Pm = Pm;
                    best_Kp = Kp;
                    best_Td = Td;
                end
            catch
                % Skip if analysis fails
                continue;
            end
        end
    end
    
    % Adapt to plant characteristics
    if plantInfo.hasRHPZeros
        % For non-minimum phase plants, reduce derivative action
        best_Td = best_Td * 0.7;
    end
    
    if plantInfo.isHighOrder
        % For high-order plants, add more filtering
        epsilon = epsilon * 1.5;
    end
    
    Kp = best_Kp;
    Kd = Kp * best_Td;
end
   
%Bis hier hat Claude

function [Kp, Ki, Kd] = findOptimalPID(G, M_s, epsilon, plantInfo)
    % Find optimal PID parameters under M_s constraint
    % Multi-stage approach for PID design
    
    % 1. Find good PI parameters
    [Kp_pi, Ki_pi] = findOptimalPI(G, M_s * 1.1, plantInfo);
    
    % 2. Add derivative action starting from PI solution
    Td_values = logspace(-2, 0.5, 10);
    
    best_Kp = Kp_pi;
    best_Ki = Ki_pi;
    best_Td = 0;
    best_score = -Inf;
    
    for j = 1:length(Td_values)
        Td = Td_values(j);
        Kd = Kp_pi * Td;
        
        K_test = tf([Kd, Kp_pi, Ki_pi], [epsilon*Td, 1, 0]);
        L = G * K_test;
        
        try
            % Calculate closed-loop response
            S = feedback(1, L);
            [peakgain, ~, ~] = getPeakGain(S);
            [~, Pm] = margin(L);
            
            if peakgain <= M_s
                % Score based on phase margin and integral gain
                score = Ki_pi * sqrt(Kd) * (Pm/100);
                
                if score > best_score
                    best_score = score;
                    best_Td = Td;
                end
            end
        catch
            % Skip if analysis fails
            continue;
        end
    end
    
    % 3. Fine-tune the PID controller
    best_Kp = Kp_pi;
    best_Ki = Ki_pi;
    best_Kd = best_Kp * best_Td;
    
    % Adjust parameters based on plant characteristics
    if plantInfo.hasIntegrator
        % For plants with integrator, reduce integral action
        best_Ki = best_Ki * 0.7;
    end
    
    if plantInfo.hasRHPZeros
        % For non-minimum phase plants, reduce derivative action
        best_Kd = best_Kd * 0.7;
    end
    
    if plantInfo.isHighOrder
        % For high-order plants, add more filtering
        epsilon = epsilon * 1.5;
    end
    
    Kp = best_Kp;
    Ki = best_Ki;
    Kd = best_Kd;
end

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

% Support function for plant information string
function infoStr = getPlantInfoString(plantInfo)
    % GETPLANTINFOSTRING Creates a formatted string with plant characteristics
    
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

% Support function for bandwidth estimation
function w_bandwidth = getBandwidth(G, plantInfo)
    % GETBANDWIDTH Estimate bandwidth of plant
    
    try
        % Try standard bandwidth calculation for stable systems
        if ~plantInfo.isUnstable
            w_bandwidth = bandwidth(G);
            if ~isnan(w_bandwidth)
                return;
            end
        end
    catch
        % Continue to alternative methods if bandwidth calculation fails
    end
    
    % Alternative methods for bandwidth estimation
    try
        % Method 1: Use frequency response
        w = logspace(-3, 3, 500);
        [mag, ~] = bode(G, w);
        mag = squeeze(mag);
        
        % Find 0.7 (-3dB) crossing point
        try
            dc_gain = dcgain(G);
        catch
            % If DC gain calculation fails, use low-frequency gain
            dc_gain = mag(1);
        end
        
        threshold = 0.7 * abs(dc_gain);
        cross_idx = find(mag < threshold, 1);
        
        if ~isempty(cross_idx) && cross_idx > 1
            w_bandwidth = w(cross_idx);
            return;
        end
    catch
        % Continue to next method
    end
    
    % Method 2: Use pole/zero information
    try
        p = plantInfo.poles;
        
        if ~isempty(p)
            % Use dominant (slowest) pole for stable systems 
            % or fastest unstable pole for unstable systems
            if plantInfo.isUnstable
                unstable_poles = p(real(p) > 0);
                if ~isempty(unstable_poles)
                    [~, idx] = max(real(unstable_poles));
                    w_bandwidth = 5 * abs(unstable_poles(idx));  % Higher bandwidth needed for unstable pole
                    return;
                end
            else
                % For stable systems, use dominant poles
                stable_poles = p(real(p) < 0);
                if ~isempty(stable_poles)
                    [~, idx] = min(abs(real(stable_poles)));
                    w_bandwidth = 2 * abs(real(stable_poles(idx)));  % Rule of thumb
                    return;
                end
            end
        end
    catch
        % Continue to next method
    end
    
    % Method 3: Use FOPDT parameters if available
    if ~isnan(plantInfo.FOPDT.T) && plantInfo.FOPDT.T > 0
        w_bandwidth = 2.5 / plantInfo.FOPDT.T;  % Rule of thumb
        return;
    end
    
    % Default fallback value
    w_bandwidth = 1.0;
end

% Hilfsfunktion zur Berechnung der maximalen Verstärkung
function [peakgain, wpeak, w] = getPeakGain(sys)
    % Berechnet den maximalen Amplitudengang eines Systems
    w = logspace(-3, 3, 1000);
    [mag, ~] = bode(sys, w);
    mag = squeeze(mag);
    [peakgain, idx] = max(mag);
    wpeak = w(idx);
end

% EVALUATECONTROLLER Evaluates a controller based on specified criteria
function score = evaluateController(K, G, goal, desired_pm, desired_os, desired_ts, desired_bw)
    % Initialize score
    score = 0;
    
    % Ensure goal is in correct format
    if strcmpi(goal, 'tracking')
        goal = 'Tracking';
    elseif strcmpi(goal, 'disturbance rejection')
        goal = 'Disturbance Rejection';
    elseif strcmpi(goal, 'robustness')
        goal = 'Robustness';
    end
    
    % Compute closed-loop transfer functions
    try
        T = feedback(G*K, 1);      % Complementary sensitivity function
        S = feedback(1, G*K);      % Sensitivity function
        L = G*K;                   % Open-loop transfer function
        
        % Test if closed-loop system is stable
        if any(real(pole(T)) > 0)
            disp('Controller results in unstable closed-loop system.');
            score = -100;  % Unstable system gets a heavily negative score
            return;
        end
    catch ME
        disp(['Error computing closed-loop transfer functions: ', ME.message]);
        score = -50;  % Error indicates likely problems with the controller
        return;
    end
    
    % 1. Stability metrics (30 points maximum)
    try
        % Compute gain and phase margins
        [gm, pm, wgm, wpm] = margin(L);
        
        % Convert gain margin to dB
        gm_dB = 20*log10(gm);
        
        % Phase margin scoring (0-15 points)
        pm_score = 15 * (1 - min(1, abs(pm - desired_pm) / 45));
        
        % Gain margin scoring (0-15 points)
        gm_score = 15 * (1 - min(1, abs(gm_dB - 10) / 10));
        
        % Display stability metrics
        disp(['Phase Margin: ', num2str(pm), '° (desired: ', num2str(desired_pm), '°)']);
        disp(['Gain Margin: ', num2str(gm_dB), ' dB']);
        
        % Add to total score
        score = score + pm_score + gm_score;
    catch ME
        disp(['Warning: Could not compute stability margins: ', ME.message]);
        % Apply penalty for not being able to compute stability margins
        score = score - 15;
    end
    
    % 2. Time-domain performance (40 points maximum)
    try
        % Compute step response
        [y, t] = step(T);
        
        % Extract step response metrics
        info = stepinfo(y, t);
        
        % Calculate key metrics
        actual_overshoot = info.Overshoot;
        actual_settling = info.SettlingTime;
        actual_rise = info.RiseTime;
        
        % Display time domain metrics
        disp(['Overshoot: ', num2str(actual_overshoot), '% (desired: ', num2str(desired_os), '%)']);
        disp(['Settling Time: ', num2str(actual_settling), ' s (desired: ', num2str(desired_ts), ' s)']);
        disp(['Rise Time: ', num2str(actual_rise), ' s']);
        
        % Score overshoot (0-15 points)
        % Better score for being close to desired overshoot
        os_score = 15 * (1 - min(1, abs(actual_overshoot - desired_os) / 50));
        
        % Score settling time (0-15 points)
        % Better score for faster settling
        ts_score = 15 * (1 - min(1, abs(actual_settling - desired_ts) / (2*desired_ts)));
        
        % Score rise time (0-10 points)
        % Better score for faster rise time relative to settling time
        rt_score = 10 * (1 - min(1, actual_rise / desired_ts));
        
        % Add to total score
        score = score + os_score + ts_score + rt_score;
    catch ME
        disp(['Warning: Could not compute time-domain performance metrics: ', ME.message]);
        % Apply penalty for not being able to compute time domain metrics
        score = score - 20;
    end
    
    % 3. Frequency-domain performance (20 points maximum)
    try
        % Calculate bandwidth
        bw = bandwidth(T);
        
        % Display bandwidth
        disp(['Bandwidth: ', num2str(bw), ' rad/s (desired: ', num2str(desired_bw), ' rad/s)']);
        
        % Score bandwidth (0-20 points)
        % Better score for being close to desired bandwidth
        bw_score = 20 * (1 - min(1, abs(bw - desired_bw) / desired_bw));
        
        % Add to total score
        score = score + bw_score;
    catch ME
        disp(['Warning: Could not compute bandwidth: ', ME.message]);
        % Apply penalty for not being able to compute bandwidth
        score = score - 10;
    end
    
    % 4. Goal-specific metrics (10 points maximum)
    try
        switch goal
            case 'Tracking'
                % For tracking, evaluate the reference tracking performance
                % Calculate integral of error for a step reference
                t_sim = linspace(0, 3*desired_ts, 1000);
                r = ones(size(t_sim));
                [y_r, ~] = lsim(T, r, t_sim);
                err_integral = trapz(t_sim, abs(r - y_r));
                
                % Score tracking (0-10 points)
                tracking_score = 10 * exp(-err_integral / 5);
                
                % Add to total score
                score = score + tracking_score;
                
                disp(['Tracking Performance Score: ', num2str(tracking_score), '/10']);
                
            case 'Disturbance Rejection'
                % For disturbance rejection, evaluate the disturbance response
                % Calculate the integral of output for a step disturbance
                t_sim = linspace(0, 3*desired_ts, 1000);
                d = ones(size(t_sim));
                [y_d, ~] = lsim(S, d, t_sim);  % S is sensitivity function
                dist_integral = trapz(t_sim, abs(y_d));
                
                % Score disturbance rejection (0-10 points)
                dist_score = 10 * exp(-dist_integral / 5);
                
                % Add to total score
                score = score + dist_score;
                
                disp(['Disturbance Rejection Score: ', num2str(dist_score), '/10']);
                
            case 'Robustness'
                % For robustness, evaluate the sensitivity peak
                try
                    [Ms, w_ms] = getPeakGain(S);
                    
                    % Good sensitivity peak should be below 2 (6 dB)
                    robust_score = 10 * exp(-(Ms-1) / 1.5);
                    
                    % Add to total score
                    score = score + robust_score;
                    
                    disp(['Sensitivity Peak (Ms): ', num2str(Ms)]);
                    disp(['Robustness Score: ', num2str(robust_score), '/10']);
                catch
                    disp('Could not compute sensitivity peak.');
                    score = score - 5;
                end
        end
    catch ME
        disp(['Warning: Could not compute goal-specific metrics: ', ME.message]);
        % Apply penalty for not being able to compute goal-specific metrics
        score = score - 5;
    end
    
    % Normalize score to be within 0-100 range
    score = max(0, min(100, score));
    
    % Display final score
    disp(['Final Controller Score: ', num2str(score), '/100']);
end
