function score = evaluateController(K, G, goal, desired_pm, desired_os, desired_ts, desired_bw, plantInfo)
    % EVALUATECONTROLLER Evaluates a controller based on specified criteria with
    % enhanced assessment for difficult plants
    %
    % Inputs:
    %   K         - Controller transfer function
    %   G         - Plant transfer function
    %   goal      - Optimization goal: 'Tracking', 'Disturbance Rejection', 'Robustness'
    %   desired_pm - Desired phase margin in degrees
    %   desired_os - Desired overshoot percentage
    %   desired_ts - Desired settling time
    %   desired_bw - Desired bandwidth
    %   plantInfo - (Optional) Plant information structure from analyzePlant
    %
    % Output:
    %   score - Controller performance score (0-100)
    
    % Create empty plantInfo if not provided
    if nargin < 8 || isempty(plantInfo)
        try
            plantInfo = analyzePlant(G);
        catch
            % Create minimal plantInfo
            plantInfo = struct();
            plantInfo.isUnstable = any(real(pole(G)) > 0);
            plantInfo.hasRHPZeros = false;
            try
                z = zero(G);
                plantInfo.hasRHPZeros = any(real(z) > 0);
            catch
                % No zeros or couldn't compute
            end
            plantInfo.isHighOrder = length(pole(G)) > 3;
            plantInfo.hasIntegrator = any(abs(pole(G)) < 1e-6);
        end
    end
    
    % Initialize score
    score = 50;  % Start from neutral position
    
    % Ensure goal is in correct format
    if strcmpi(goal, 'tracking')
        goal = 'Tracking';
    elseif strcmpi(goal, 'disturbance rejection') || strcmpi(goal, 'disturbance')
        goal = 'Disturbance Rejection';
    elseif strcmpi(goal, 'robustness')
        goal = 'Robustness';
    end
    
    % Adaptive weighting factors based on plant characteristics
    weights = getAdaptiveWeights(plantInfo, goal);
    
    % Compute closed-loop transfer functions
    try
        % Standard transfer functions
        L = G*K;                    % Open-loop transfer function
        T = feedback(L, 1);         % Complementary sensitivity
        S = feedback(1, L);         % Sensitivity
        CS = minreal(K*S);          % Control sensitivity (noise sensitivity)
        GS = minreal(G*S);          % Load disturbance sensitivity
        
        % Test if closed-loop system is stable
        if any(real(pole(T)) > 0)
            disp('Controller results in unstable closed-loop system.');
            score = -100;  % Unstable system gets a heavily negative score
            return;
        end
    catch ME
        disp(['Error computing closed-loop transfer functions: ', ME.message]);
        score = 0;  % Reset score for proper evaluation
        return;
    end
    
    % 1. Stability metrics
    stabScore = evaluateStabilityMetrics(L, K, G, desired_pm, weights.stability);
    
    % 2. Time-domain performance 
    timeScore = evaluateTimeDomainPerformance(T, desired_os, desired_ts, weights.time);
    
    % 3. Frequency-domain performance
    freqScore = evaluateFrequencyDomainPerformance(T, S, CS, desired_bw, weights.frequency);
    
    % 4. Goal-specific metrics
    goalScore = evaluateGoalSpecificMetrics(T, S, CS, GS, goal, desired_ts, weights.goal);
    
    % 5. Controller complexity and practical aspects
    complexityScore = evaluateControllerComplexity(K, G, weights.complexity);
    
    % Calculate final score
    score = stabScore + timeScore + freqScore + goalScore + complexityScore;
    
    % Normalize to 0-100 range
    score = min(max(score, 0), 100);
    
    % Display comprehensive results if not suppressed
    displayResults(stabScore, timeScore, freqScore, goalScore, complexityScore, score);
end

function weights = getAdaptiveWeights(plantInfo, goal)
    % Determine weights for different evaluation components based on plant difficulty
    weights = struct();
    
    % Base weights
    weights.stability = 30;
    weights.time = 25;
    weights.frequency = 20;
    weights.goal = 15;
    weights.complexity = 10;
    
    % Adjust for difficult plants
    if plantInfo.isUnstable
        % For unstable plants, stability is more important
        weights.stability = 40;
        weights.time = 20;
        weights.frequency = 15;
        weights.goal = 15;
        weights.complexity = 10;
    end
    
    if plantInfo.hasRHPZeros
        % For non-minimum phase plants, tradeoffs are more important
        weights.frequency = weights.frequency + 5;
        weights.time = weights.time - 5;
    end
    
    if plantInfo.isHighOrder
        % For high-order plants, controller complexity matters more
        weights.complexity = weights.complexity + 5;
        weights.frequency = weights.frequency - 5;
    end
    
    % Adjust for specific goals
    if strcmp(goal, 'Robustness')
        weights.stability = weights.stability + 5;
        weights.goal = weights.goal + 5;
        weights.time = weights.time - 5;
        weights.frequency = weights.frequency - 5;
    elseif strcmp(goal, 'Tracking')
        weights.time = weights.time + 5;
        weights.stability = weights.stability - 5;
    elseif strcmp(goal, 'Disturbance Rejection')
        weights.goal = weights.goal + 5;
        weights.complexity = weights.complexity - 5;
    end
    
    % Normalize weights to sum to 100
    total = weights.stability + weights.time + weights.frequency + weights.goal + weights.complexity;
    factor = 100 / total;
    
    weights.stability = weights.stability * factor;
    weights.time = weights.time * factor;
    weights.frequency = weights.frequency * factor;
    weights.goal = weights.goal * factor;
    weights.complexity = weights.complexity * factor;
end

function score = evaluateStabilityMetrics(L, K, G, desired_pm, maxScore)
    % Evaluate stability characteristics of the control system
    score = 0;
    try
        % Compute gain and phase margins
        [gm, pm, wgm, wpm] = margin(L);
        
        % Convert gain margin to dB
        gm_dB = 20*log10(gm);
        
        % Phase margin scoring (0-40% of maxScore)
        pm_target = desired_pm;
        if isnan(pm_target) || pm_target < 30
            pm_target = 45; % Default if not specified or too low
        end
        
        pm_score = 0.4 * maxScore * (1 - min(1, abs(pm - pm_target) / 45));
        
        % Gain margin scoring (0-30% of maxScore)
        gm_target = 10; % 10 dB is a good target
        gm_score = 0.3 * maxScore * (1 - min(1, abs(gm_dB - gm_target) / 10));
        
        % Sensitivity peak scoring (0-30% of maxScore)
        try
            S = feedback(1, L);
            [Ms, wMs] = getPeakGain(S);
            
            % Sensitivity peak should be less than 2 (6dB) for good robustness
            Ms_score = 0.3 * maxScore * exp(-(Ms-1) / 1.5);
            
            disp(['Sensitivity Peak (Ms): ', num2str(Ms), ' at ', num2str(wMs), ' rad/s']);
        catch
            % If sensitivity peak calculation fails
            Ms_score = 0;
            disp('Could not compute sensitivity peak.');
        end
        
        % Display stability metrics
        disp(['Phase Margin: ', num2str(pm), '° at ', num2str(wpm), ' rad/s']);
        disp(['Gain Margin: ', num2str(gm_dB), ' dB at ', num2str(wgm), ' rad/s']);
        
        % Add internal stability check for unstable plants
        if any(real(pole(G)) > 0) || any(real(pole(K)) > 0)
            % Check for unstable pole-zero cancellations
            try
                % Check internal stability using additional test
                T = feedback(L, 1);
                KS = K * feedback(1, L);
                GS = G * feedback(1, L);
                
                internal_stable = all(real(pole(T)) < 0) && all(real(pole(KS)) < 0) && all(real(pole(GS)) < 0);
                
                if internal_stable
                    % Bonus for properly stabilizing an unstable plant
                    internal_score = 0.1 * maxScore;
                    disp('Internal stability verified for unstable plant - Good!');
                else
                    % Major penalty for internal stability issues
                    internal_score = -0.5 * maxScore;
                    disp('WARNING: Internal stability issues detected!');
                end
            catch
                internal_score = 0;
                disp('Could not verify internal stability.');
            end
        else
            internal_score = 0;
        end
        
        % Compute total stability score
        score = pm_score + gm_score + Ms_score + internal_score;
        
        % Ensure score is within bounds for this category
        score = min(max(score, 0), maxScore);
        
    catch ME
        disp(['Warning: Could not compute stability margins: ', ME.message]);
        score = 0.1 * maxScore; % Small base score if stable but margins can't be computed
    end
end

function score = evaluateTimeDomainPerformance(T, desired_os, desired_ts, maxScore)
    % Evaluate time-domain performance
    score = 0;
    try
        % Compute step response with sufficient time for settling
        sim_time = max(20, 3*desired_ts);
        [y, t] = step(T, linspace(0, sim_time, 1000));
        
        % Extract step response metrics
        info = stepinfo(y, t);
        
        % Calculate key metrics
        actual_overshoot = info.Overshoot;
        actual_settling = info.SettlingTime;
        actual_rise = info.RiseTime;
        
        % Check final value (steady-state error)
        if t(end) > 5*desired_ts
            final_value = y(end);
            ss_error = abs(1 - final_value);
            
            % Score steady-state error (0-25% of maxScore)
            if ss_error < 0.01
                ss_score = 0.25 * maxScore;
            elseif ss_error < 0.05
                ss_score = 0.15 * maxScore;
            elseif ss_error < 0.1
                ss_score = 0.05 * maxScore;
            else
                ss_score = 0;
            end
        else
            ss_score = 0.1 * maxScore; % Default if simulation time not sufficient
            ss_error = NaN;
        end
        
        % Display time domain metrics
        disp(['Overshoot: ', num2str(actual_overshoot), '% (desired: ', num2str(desired_os), '%)']);
        disp(['Settling Time: ', num2str(actual_settling), ' s (desired: ', num2str(desired_ts), ' s)']);
        disp(['Rise Time: ', num2str(actual_rise), ' s']);
        disp(['Steady-State Error: ', num2str(ss_error * 100), '%']);
        
        % Score overshoot (0-30% of maxScore)
        % Better score for being close to desired overshoot
        if isnan(actual_overshoot)
            os_score = 0.1 * maxScore; % Default if not computable
        else
            os_target = desired_os;
            if isnan(os_target)
                os_target = 10; % Default if not specified
            end
            
            % Zero overshoot is acceptable even if some overshoot is desired
            if actual_overshoot == 0 && os_target < 20
                os_score = 0.2 * maxScore;
            else
                os_score = 0.3 * maxScore * (1 - min(1, abs(actual_overshoot - os_target) / 30));
            end
        end
        
        % Score settling time (0-30% of maxScore)
        % Better score for faster settling
        if isnan(actual_settling)
            ts_score = 0.1 * maxScore; % Default if not computable
        else
            ts_score = 0.3 * maxScore * (1 - min(1, abs(actual_settling - desired_ts) / (2*desired_ts)));
        end
        
        % Score rise time (0-15% of maxScore)
        % Better score for faster rise time relative to settling time
        if isnan(actual_rise)
            rt_score = 0.05 * maxScore; % Default if not computable
        else
            rt_score = 0.15 * maxScore * (1 - min(1, actual_rise / desired_ts));
        end
        
        % Compute total time-domain score
        score = os_score + ts_score + rt_score + ss_score;
        
        % Ensure score is within bounds for this category
        score = min(max(score, 0), maxScore);
        
    catch ME
        disp(['Warning: Could not compute time-domain performance metrics: ', ME.message]);
        score = 0.1 * maxScore; % Small score if time-domain analysis fails
    end
end

function score = evaluateFrequencyDomainPerformance(T, S, CS, desired_bw, maxScore)
    % Evaluate frequency-domain performance
    score = 0;
    try
        % Calculate bandwidth
        bw = bandwidth(T);
        
        % Display bandwidth
        disp(['Bandwidth: ', num2str(bw), ' rad/s (desired: ', num2str(desired_bw), ' rad/s)']);
        
        % Score bandwidth (0-40% of maxScore)
        % Better score for being close to desired bandwidth
        bw_score = 0.4 * maxScore * (1 - min(1, abs(bw - desired_bw) / desired_bw));
        
        % Evaluate high-frequency rolloff (0-30% of maxScore)
        try
            w = logspace(log10(bw), log10(bw) + 2, 100);
            [mag_t, ~] = bode(T, w);
            mag_t = squeeze(mag_t);
            
            % Calculate rolloff rate (dB/decade)
            idx1 = 1;
            idx2 = length(w);
            rolloff = -20 * log10(mag_t(idx2)/mag_t(idx1)) / (log10(w(idx2)/w(idx1)));
            
            % Good rolloff is at least 20 dB/decade, excellent is 40+ dB/decade
            if rolloff >= 40
                rolloff_score = 0.3 * maxScore;
            elseif rolloff >= 20
                rolloff_score = 0.2 * maxScore;
            elseif rolloff >= 10
                rolloff_score = 0.1 * maxScore;
            else
                rolloff_score = 0;
            end
            
            disp(['High-frequency rolloff: ', num2str(rolloff), ' dB/decade']);
        catch
            rolloff_score = 0.1 * maxScore; % Default if rolloff can't be computed
        end
        
        % Evaluate noise sensitivity (0-30% of maxScore)
        try
            [mag_cs, ~] = bode(CS, w);
            mag_cs = squeeze(mag_cs);
            max_cs = max(mag_cs);
            
            % Control sensitivity should be low at high frequencies to reduce noise amplification
            if max_cs < 5
                noise_score = 0.3 * maxScore;
            elseif max_cs < 10
                noise_score = 0.2 * maxScore;
            elseif max_cs < 20
                noise_score = 0.1 * maxScore;
            else
                noise_score = 0;
            end
            
            disp(['Noise Sensitivity Peak: ', num2str(max_cs)]);
        catch
            noise_score = 0.1 * maxScore; % Default if noise sensitivity can't be computed
        end
        
        % Compute total frequency-domain score
        score = bw_score + rolloff_score + noise_score;
        
        % Ensure score is within bounds for this category
        score = min(max(score, 0), maxScore);
        
    catch ME
        disp(['Warning: Could not compute frequency-domain performance metrics: ', ME.message]);
        score = 0.1 * maxScore; % Small score if frequency-domain analysis fails
    end
end

function score = evaluateGoalSpecificMetrics(T, S, CS, GS, goal, desired_ts, maxScore)
    % Evaluate performance metrics specific to the control goal
    score = 0;
    try
        switch goal
            case 'Tracking'
                % For tracking, evaluate the reference tracking performance
                % Calculate integral of error for a step reference
                t_sim = linspace(0, 5*desired_ts, 1000);
                r = ones(size(t_sim));
                try
                    [y_r, ~] = lsim(T, r, t_sim);
                    err_integral = trapz(t_sim, abs(r - y_r));
                    
                    % Score tracking (0-70% of maxScore)
                    tracking_score = 0.7 * maxScore * exp(-err_integral / 3);
                    
                    disp(['Tracking IAE: ', num2str(err_integral)]);
                catch
                    tracking_score = 0.3 * maxScore; % Default if simulation fails
                end
                
                % Evaluate input usage (0-30% of maxScore)
                try
                    Ksys = CS; % CS = KS, control sensitivity function
                    [u, ~] = lsim(Ksys, r, t_sim);
                    
                    % Calculate control effort metrics
                    max_u = max(abs(u));
                    u_integral = trapz(t_sim, abs(u));
                    
                    % Reasonable control effort is desired
                    if max_u < 10 && u_integral < 20*desired_ts
                        control_score = 0.3 * maxScore;
                    elseif max_u < 20 && u_integral < 40*desired_ts
                        control_score = 0.2 * maxScore;
                    elseif max_u < 50
                        control_score = 0.1 * maxScore;
                    else
                        control_score = 0;
                    end
                    
                    disp(['Maximum Control Effort: ', num2str(max_u)]);
                    disp(['Control Effort Integral: ', num2str(u_integral)]);
                catch
                    control_score = 0.1 * maxScore; % Default if simulation fails
                end
                
                % Total score for tracking goal
                score = tracking_score + control_score;
                
            case 'Disturbance Rejection'
                % For disturbance rejection, evaluate the disturbance response
                % Calculate the integral of output for a step disturbance
                t_sim = linspace(0, 5*desired_ts, 1000);
                d = ones(size(t_sim));
                
                try
                    [y_d, ~] = lsim(S, d, t_sim);  % S is sensitivity function for output disturbance
                    dist_integral = trapz(t_sim, abs(y_d));
                    
                    % Score disturbance rejection (0-60% of maxScore)
                    dist_score = 0.6 * maxScore * exp(-dist_integral / 3);
                    
                    disp(['Disturbance Rejection IAE: ', num2str(dist_integral)]);
                catch
                    dist_score = 0.3 * maxScore; % Default if simulation fails
                end
                
                % Evaluate load disturbance rejection
                try
                    [y_load, ~] = lsim(GS, d, t_sim);  % GS for input/load disturbance
                    load_integral = trapz(t_sim, abs(y_load));
                    
                    % Score load disturbance rejection (0-40% of maxScore)
                    load_score = 0.4 * maxScore * exp(-load_integral / 3);
                    
                    disp(['Load Disturbance Rejection IAE: ', num2str(load_integral)]);
                catch
                    load_score = 0.2 * maxScore; % Default if simulation fails
                end
                
                % Total score for disturbance rejection goal
                score = dist_score + load_score;
                
            case 'Robustness'
                % For robustness, evaluate multiple sensitivity functions
                
                % Sensitivity peak (0-40% of maxScore)
                try
                    [Ms, w_ms] = getPeakGain(S);
                    
                    % Good sensitivity peak should be below 2 (6 dB)
                    S_score = 0.4 * maxScore * exp(-(Ms-1) / 1.2);
                    
                    disp(['Sensitivity Peak (Ms): ', num2str(Ms), ' at ', num2str(w_ms), ' rad/s']);
                catch
                    S_score = 0.2 * maxScore; % Default if sensitivity peak can't be computed
                end
                
                % Complementary sensitivity peak (0-30% of maxScore)
                try
                    [Mt, w_mt] = getPeakGain(T);
                    
                    % Good complementary sensitivity peak should be below 1.5 (3.5 dB)
                    T_score = 0.3 * maxScore * exp(-(Mt-1) / 1);
                    
                    disp(['Complementary Sensitivity Peak (Mt): ', num2str(Mt), ' at ', num2str(w_mt), ' rad/s']);
                catch
                    T_score = 0.15 * maxScore; % Default if complementary sensitivity peak can't be computed
                end
                
                % Gain variation tolerance (0-30% of maxScore)
                try
                    % Evaluate complementary sensitivity at low frequencies
                    w_low = logspace(-3, -1, 50);
                    [mag_t_low, ~] = bode(T, w_low);
                    mag_t_low = squeeze(mag_t_low);
                    
                    % Evaluate sensitivity at high frequencies
                    w_high = logspace(1, 3, 50);
                    [mag_s_high, ~] = bode(S, w_high);
                    mag_s_high = squeeze(mag_s_high);
                    
                    % Good robustness to gain variations
                    if max(mag_t_low) < 1.1 && max(mag_s_high) < 1.1
                        gain_score = 0.3 * maxScore;
                    elseif max(mag_t_low) < 1.3 && max(mag_s_high) < 1.3
                        gain_score = 0.2 * maxScore;
                    elseif max(mag_t_low) < 1.5 && max(mag_s_high) < 1.5
                        gain_score = 0.1 * maxScore;
                    else
                        gain_score = 0;
                    end
                catch
                    gain_score = 0.1 * maxScore; % Default if gain variation analysis fails
                end
                
                % Total score for robustness goal
                score = S_score + T_score + gain_score;
                
            otherwise
                % Unknown goal, use general performance metrics
                score = 0.5 * maxScore;
        end
        
        % Ensure score is within bounds for this category
        score = min(max(score, 0), maxScore);
        
    catch ME
        disp(['Warning: Could not compute goal-specific metrics: ', ME.message]);
        score = 0.1 * maxScore; % Small score if goal-specific analysis fails
    end
end

function score = evaluateControllerComplexity(K, G, maxScore)
    % Evaluate controller complexity and practical aspects
    score = 0;
    try
        % Get controller order
        [num, den] = tfdata(K, 'v');
        ctrl_order = length(den) - 1;
        
        % Get plant order
        [num_g, den_g] = tfdata(G, 'v');
        plant_order = length(den_g) - 1;
        
        % Calculate relative controller complexity
        relative_complexity = ctrl_order / max(1, plant_order);
        
        % Determine controller structure
        ctrl_type = determineControllerType(K);
        
        % Display controller information
        disp(['Controller Type: ', ctrl_type]);
        disp(['Controller Order: ', num2str(ctrl_order)]);
        disp(['Relative Complexity: ', num2str(relative_complexity)]);
        
        % Score controller structure (0-40% of maxScore)
        switch ctrl_type
            case 'P'
                structure_score = 0.4 * maxScore;
            case 'PI'
                structure_score = 0.35 * maxScore;
            case 'PD'
                structure_score = 0.3 * maxScore;
            case 'PID'
                structure_score = 0.25 * maxScore;
            otherwise
                % More complex controllers get lower scores
                structure_score = 0.2 * maxScore * (1 - min(0.8, relative_complexity/5));
        end
        
        % Score controller coefficient magnitudes (0-30% of maxScore)
        % Check for reasonable coefficient magnitudes
        if max(abs([num, den])) < 1e6 && min(abs([num(num~=0), den(den~=0)])) > 1e-6
            coeff_score = 0.3 * maxScore;
        elseif max(abs([num, den])) < 1e8 && min(abs([num(num~=0), den(den~=0)])) > 1e-8
            coeff_score = 0.2 * maxScore;
        else
            coeff_score = 0.1 * maxScore;
        end
        
        % Score pole-zero properties (0-30% of maxScore)
        try
            z = zero(K);
            p = pole(K);
            
            % Check for controller stability
            if all(real(p) < 0)
                % Stable controller is good
                pz_score = 0.15 * maxScore;
                
                % Check for reasonable pole-zero spacing
                min_spacing = Inf;
                if ~isempty(z) && ~isempty(p)
                    for i = 1:length(p)
                        for j = 1:length(z)
                            spacing = abs(p(i) - z(j));
                            min_spacing = min(min_spacing, spacing);
                        end
                    end
                    
                    % Good pole-zero spacing avoids nearly cancelling dynamics
                    if min_spacing > 0.1
                        pz_score = pz_score + 0.15 * maxScore;
                    elseif min_spacing > 0.01
                        pz_score = pz_score + 0.05 * maxScore;
                    end
                else
                    % If no zeros or poles, full score
                    pz_score = 0.3 * maxScore;
                end
            else
                % Unstable controller is generally bad
                pz_score = 0;
            end
        catch
            % Default if pole-zero analysis fails
            pz_score = 0.1 * maxScore;
        end
        
        % Compute total complexity score
        score = structure_score + coeff_score + pz_score;
        
        % Ensure score is within bounds for this category
        score = min(max(score, 0), maxScore);
        
    catch ME
        disp(['Warning: Could not evaluate controller complexity: ', ME.message]);
        score = 0.1 * maxScore; % Small score if complexity analysis fails
    end
end

function displayResults(stabScore, timeScore, freqScore, goalScore, complexityScore, totalScore)
    % Display comprehensive results summary
    disp('=== CONTROLLER EVALUATION SUMMARY ===');
    disp(['Stability Metrics Score:    ', num2str(stabScore, '%.1f')]);
    disp(['Time-Domain Performance:    ', num2str(timeScore, '%.1f')]);
    disp(['Frequency-Domain Performance: ', num2str(freqScore, '%.1f')]);
    disp(['Goal-Specific Metrics:      ', num2str(goalScore, '%.1f')]);
    disp(['Controller Complexity:      ', num2str(complexityScore, '%.1f')]);
    disp('-----------------------------------');
    disp(['TOTAL CONTROLLER SCORE:     ', num2str(totalScore, '%.1f'), '/100']);
    
    % Interpret score
    if totalScore >= 90
        disp('RATING: EXCELLENT - Controller exceeds all design requirements');
    elseif totalScore >= 80
        disp('RATING: VERY GOOD - Controller meets all key requirements with excellent performance');
    elseif totalScore >= 70
        disp('RATING: GOOD - Controller meets requirements with good performance');
    elseif totalScore >= 60
        disp('RATING: SATISFACTORY - Controller meets basic requirements');
    elseif totalScore >= 50
        disp('RATING: ADEQUATE - Controller is functional but has limitations');
    elseif totalScore >= 30
        disp('RATING: POOR - Controller has significant performance issues');
    else
        disp('RATING: UNSATISFACTORY - Controller fails to meet basic requirements');
    end
end

function [peakgain, wpeak, w] = getPeakGain(sys)
    % Get peak gain of a system and the frequency where it occurs
    w = logspace(-3, 3, 500);
    [mag, ~] = bode(sys, w);
    mag = squeeze(mag);
    
    % Find peak
    [peakgain, idx] = max(mag);
    wpeak = w(idx);
end

function controllerType = determineControllerType(K)
    % DETERMINECONTROLLERTYPE - Determine controller type from transfer function
    
    % Get numerator and denominator
    [num, den] = tfdata(K, 'v');
    
    % Remove leading zeros
    num = num(find(abs(num) > 1e-10, 1):end);
    den = den(find(abs(den) > 1e-10, 1):end);
    
    % Check for different controller types
    has_integrator = any(abs(den) < 1e-10);
    
    % Determine type based on transfer function structure
    if length(num) == 1 && length(den) == 1
        % Simple proportional controller
        controllerType = 'P';
    elseif has_integrator && length(num) <= 2 && length(den) <= 2
        % PI controller
        controllerType = 'PI';
    elseif ~has_integrator && length(num) >= 2 && length(den) >= 2 && length(num) <= 3 && length(den) <= 3
        % PD controller with filter
        controllerType = 'PD';
    elseif has_integrator && length(num) >= 3 && length(den) >= 2 && length(den) <= 3
        % PID controller
        controllerType = 'PID';
    elseif length(num) > 3 || length(den) > 3
        % Higher-order controller
        controllerType = 'Higher-Order';
    else
        % Default to unknown type
        controllerType = 'Custom';
    end
end