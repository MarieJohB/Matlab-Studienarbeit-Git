function [K, details] = designCompensationController(G, structure, options, plantInfo)
% DESIGNCOMPENSATIONCONTROLLER Design a controller that directly compensates for problematic plant dynamics
% Enhanced version with better handling of highly unstable, non-minimum phase, and high-order plants
%
% Inputs:
%   G         - Plant model (transfer function)
%   structure - Controller structure ('P', 'PI', 'PD', 'PID')
%   options   - Design options structure
%   plantInfo - Plant analysis information
%
% Outputs:
%   K        - Designed controller
%   details  - Text description of the design process

    % Initialize details
    details = 'COMPENSATION CONTROLLER DESIGN\n';
    details = [details, '----------------------------\n'];
    details = [details, sprintf('Plant Analysis: %s\n\n', getPlantInfoString(plantInfo))];
    
    % Extract plant poles and zeros
    p = plantInfo.poles;
    z = plantInfo.zeros;
    
    % Extract options with defaults
    bandwidth = getOption(options, 'bandwidth', 1);
    damping = getOption(options, 'damping', 0.8);
    epsilon = getOption(options, 'epsilon', 0.1);
    
    % Analyze plant to detect special cases
    details = [details, '1. Analyzing Plant Dynamics:\n'];
    
    % Handle very challenging plants with special analysis
    [isVeryUnstable, isHighlyOscillatory, hasCloseRHPZeros, plantCategory] = analyzePlantDifficulty(plantInfo);
    
    % Adjust design parameters based on plant difficulty
    if isVeryUnstable
        details = [details, '   - Plant has multiple/severe unstable poles - using conservative design\n'];
        bandwidth = min(bandwidth, getUnstableBandwidthLimit(plantInfo));
        damping = max(damping, 0.9); % More damping for stability
        epsilon = max(epsilon, 0.2); % More filtering for noise rejection
    end
    
    if hasCloseRHPZeros
        details = [details, '   - Plant has RHP zeros close to imaginary axis - limiting bandwidth\n'];
        bandwidth = min(bandwidth, getRHPZeroBandwidthLimit(plantInfo));
    end
    
    if isHighlyOscillatory
        details = [details, '   - Plant has poorly damped oscillatory modes - adding damping compensation\n'];
    end
    
    details = [details, sprintf('   - Plant category: %s\n', plantCategory)];
    details = [details, sprintf('   - Adjusted bandwidth: %.3f rad/s\n', bandwidth)];
    details = [details, sprintf('   - Adjusted damping: %.3f\n', damping)];
    details = [details, sprintf('   - Adjusted epsilon: %.3f\n', epsilon)];
    
    try
        % Get plant in zero-pole-gain form
        [zeros_G, poles_G, gain_G] = zpkdata(G, 'v');
        
        % Step 2: Identify problematic dynamics to compensate
        unstable_poles = poles_G(real(poles_G) > 0);
        rhp_zeros = zeros_G(real(zeros_G) > 0);
        oscillatory_poles = poles_G(abs(imag(poles_G)) > 0.1*abs(real(poles_G)) & real(poles_G) < 0);
        poorly_damped = findPoorlyDampedPoles(poles_G);
        
        % For multi-unstable plants, sort by real part to handle most unstable first
        if length(unstable_poles) > 1
            [~, idx] = sort(real(unstable_poles), 'descend');
            unstable_poles = unstable_poles(idx);
        end

        details = [details, '\n2. Problematic Dynamics Identified:\n'];
        if ~isempty(unstable_poles)
            details = [details, sprintf('   - %d unstable pole(s)\n', length(unstable_poles))];
            for i = 1:min(length(unstable_poles), 3) % Show up to 3 poles
                if imag(unstable_poles(i)) ~= 0
                    details = [details, sprintf('     * %.3f + %.3fj\n', real(unstable_poles(i)), imag(unstable_poles(i)))];
                else
                    details = [details, sprintf('     * %.3f\n', real(unstable_poles(i)))];
                end
            end
        end
        if ~isempty(rhp_zeros)
            details = [details, sprintf('   - %d RHP zero(s)\n', length(rhp_zeros))];
            for i = 1:min(length(rhp_zeros), 3) % Show up to 3 zeros
                if imag(rhp_zeros(i)) ~= 0
                    details = [details, sprintf('     * %.3f + %.3fj\n', real(rhp_zeros(i)), imag(rhp_zeros(i)))];
                else
                    details = [details, sprintf('     * %.3f\n', real(rhp_zeros(i)))];
                end
            end
        end
        if ~isempty(oscillatory_poles)
            details = [details, sprintf('   - %d oscillatory pole(s)\n', length(oscillatory_poles)/2)];
        end
        if ~isempty(poorly_damped)
            details = [details, sprintf('   - %d poorly damped pole(s)\n', length(poorly_damped)/2)];
            for i = 1:2:min(length(poorly_damped), 5) % Show up to 3 pairs
                pole_i = poorly_damped(i);
                mag = abs(pole_i);
                current_damping = -real(pole_i) / mag;
                details = [details, sprintf('     * %.3f + %.3fj (damping ratio: %.3f)\n', real(pole_i), imag(pole_i), current_damping)];
            end
        end
        
        % Step 3: Design compensator components
        num_K = 1;
        den_K = 1;
        
        % Handle unstable poles with enhanced approach
        if ~isempty(unstable_poles)
            details = [details, '\n3. Compensating Unstable Poles:\n'];
            
            % For plants with multiple unstable poles, use a more conservative approach
            conservativeness = min(1.0, 0.5 + 0.25 * length(unstable_poles));
            details = [details, sprintf('   - Using conservativeness factor: %.2f\n', conservativeness)];
            
            for i = 1:length(unstable_poles)
                pole_i = unstable_poles(i);
                
                if imag(pole_i) ~= 0
                    % Skip conjugate pairs (we'll handle both at once)
                    if i < length(unstable_poles) && abs(pole_i - conj(unstable_poles(i+1))) < 1e-6
                        continue;
                    end
                    
                    if imag(pole_i) > 0
                        % For complex poles, create quadratic terms
                        real_part = real(pole_i);
                        imag_part = imag(pole_i);
                        
                        % Create numerator to cancel the unstable poles
                        quad_term_num = [1, -2*real_part, real_part^2 + imag_part^2];
                        
                        % Create stable replacement with same frequency but more damping
                        % Use negative real part with increased magnitude and preserve the imaginary part
                        new_real_part = -abs(real_part) * (1 + conservativeness); % Increase damping
                        new_quad_term = [1, -2*new_real_part, new_real_part^2 + imag_part^2];
                        
                        % Filter for more robust cancellation at higher frequencies
                        if bandwidth > 1 && abs(real_part) > 1
                            filter_freq = bandwidth * 2;
                            filter_term = [1/filter_freq^2, sqrt(2)/filter_freq, 1];
                            den_K = conv(den_K, filter_term);
                            details = [details, sprintf('   - Added high-frequency filtering at %.2f rad/s\n', filter_freq)];
                        end
                        
                        num_K = conv(num_K, quad_term_num);
                        den_K = conv(den_K, new_quad_term);
                        
                        details = [details, sprintf('   - Complex pole at %.3f+%.3fj replaced with %.3f+%.3fj\n', real_part, imag_part, new_real_part, imag_part)];
                    end
                else
                    % For real unstable poles
                    % Create numerator to cancel the unstable pole
                    num_K = conv(num_K, [1, -pole_i]);
                    
                    % Create a stable replacement (negative real part with increased magnitude)
                    stable_pole = -abs(pole_i) * (1 + conservativeness);
                    den_K = conv(den_K, [1, -stable_pole]);
                    
                    details = [details, sprintf('   - Real pole at %.3f replaced with %.3f\n', pole_i, stable_pole)];
                end
            end
        end
        
        % Handle poorly damped poles with enhanced approach
        if ~isempty(poorly_damped) && (isempty(unstable_poles) || length(unstable_poles) <= 1)
            details = [details, '\n4. Improving Poorly Damped Poles:\n'];
            
            % Sort by damping ratio (fix least damped first)
            [damping_values, sort_idx] = getPoleDamping(poorly_damped);
            poorly_damped = poorly_damped(sort_idx);
            
            % Process significant poorly damped poles
            poles_to_process = min(length(poorly_damped), 6); % Process at most 6 poles (3 pairs)
            processed = 0;
            
            for i = 1:length(poorly_damped)
                pole_i = poorly_damped(i);
                
                % Skip conjugate pairs (we'll handle both at once)
                if i < length(poorly_damped) && abs(pole_i - conj(poorly_damped(i+1))) < 1e-6
                    continue;
                end
                
                % Skip if we've processed enough poles already
                processed = processed + 2; % Count the pair
                if processed > poles_to_process
                    break;
                end
                
                if imag(pole_i) > 0
                    % For complex poles, create quadratic terms
                    real_part = real(pole_i);
                    imag_part = imag(pole_i);
                    magnitude = abs(pole_i);
                    
                    % Calculate current damping ratio
                    current_damping = -real_part / magnitude;
                    
                    if current_damping < 0.3  % Only improve damping if it's too low
                        % Create numerator to cancel the poorly damped poles
                        quad_term_num = [1, -2*real_part, real_part^2 + imag_part^2];
                        
                        % Create better damped replacement with same frequency
                        new_damping = min(0.7, current_damping * 2.5);  % More aggressive damping
                        new_real_part = -new_damping * magnitude;
                        new_quad_term = [1, -2*new_real_part, magnitude^2];
                        
                        num_K = conv(num_K, quad_term_num);
                        den_K = conv(den_K, new_quad_term);
                        
                        details = [details, sprintf('   - Increased damping for pole at %.3f+%.3fj from %.2f to %.2f\n', real_part, imag_part, current_damping, new_damping)];
                    end
                end
            end
        end
        
        % Handle RHP zeros with enhanced approach
        if ~isempty(rhp_zeros)
            details = [details, '\n5. Addressing RHP Zeros:\n'];
            details = [details, '   - RHP zeros cannot be cancelled directly\n'];
            
            % Sort zeros by real part (handle zeros closest to imaginary axis first)
            [~, idx] = sort(real(rhp_zeros));
            critical_zeros = rhp_zeros(idx);
            
            % Get the zero with smallest real part (closest to imaginary axis)
            critical_zero = critical_zeros(1);
            
            if imag(critical_zero) ~= 0
                % For complex RHP zeros
                real_part = real(critical_zero);
                filter_freq = real_part * 0.4;  % Conservative bandwidth limit
                details = [details, sprintf('   - Complex RHP zero at %.3f+%.3fj limits bandwidth\n', real(critical_zero), imag(critical_zero))];
            else
                % For real RHP zeros
                filter_freq = real(critical_zero) * 0.4;  % Conservative bandwidth limit
                details = [details, sprintf('   - Real RHP zero at %.3f limits bandwidth\n', real(critical_zero))];
            end
            
            % Add appropriate filtering based on the zero pattern
            if length(rhp_zeros) > 1
                % For multiple RHP zeros, use higher-order filter for more attenuation
                details = [details, sprintf('   - Multiple RHP zeros detected, using higher-order filter\n')];
                filter_order = min(4, length(rhp_zeros) + 1);
                
                for i = 1:filter_order/2
                    damping_factor = 0.7 + 0.1*i;  % Increase damping for each filter stage
                    filter_term = [1/filter_freq^2, 2*damping_factor/filter_freq, 1];
                    den_K = conv(den_K, filter_term);
                end
                
                details = [details, sprintf('   - Added %d-order low-pass filter with cutoff at %.3f rad/s\n', filter_order, filter_freq)];
            else
                % For a single RHP zero, use standard second-order filter
                filter_term = [1/filter_freq^2, sqrt(2)/filter_freq, 1];
                den_K = conv(den_K, filter_term);
                
                details = [details, sprintf('   - Added second-order low-pass filter with cutoff at %.3f rad/s\n', filter_freq)];
            end
        end
        
        % Step 4: Based on requested structure, add appropriate terms
        details = [details, '\n6. Applying Requested Controller Structure:\n'];
        
        % Calculate gain based on desired bandwidth and plant characteristics
        target_gain = calculateTargetGain(G, plantInfo, bandwidth);
        
        details = [details, sprintf('   - Target loop gain at crossover: %.4f\n', target_gain)];
        
        % Apply structure-specific controller elements with enhanced handling
        switch structure
            case 'P'
                Kp = target_gain;
                
                if isVeryUnstable
                    % For very unstable plants, reduce gain for safety
                    Kp = Kp * 0.5;
                    details = [details, '   - Reducing gain for highly unstable plant\n'];
                end
                
                num_K = num_K * Kp;
                
                details = [details, sprintf('   - P controller with Kp = %.4f\n', Kp)];
                
            case 'PI'
                Kp = target_gain;
                
                % Adjust for plant characteristics
                if plantInfo.hasIntegrator
                    % For plants with integrator, use very small Ki
                    Ki = Kp * bandwidth / 50;
                    details = [details, '   - Plant has integrator - using very small Ki\n'];
                elseif isVeryUnstable
                    % For very unstable plants, reduce integral action
                    Ki = Kp * bandwidth / 30;
                    details = [details, '   - Reducing integral action for highly unstable plant\n'];
                elseif hasCloseRHPZeros
                    % For plants with RHP zeros, reduce integral action
                    Ki = Kp * bandwidth / 20;
                    details = [details, '   - Reducing integral action due to RHP zeros\n'];
                else
                    % Standard calculation
                    Ki = Kp * bandwidth / 10;
                end
                
                % Create PI controller
                pi_num = [Kp, Ki];
                pi_den = [1, 0];
                
                % Combine with numerator and denominator
                num_K = conv(num_K, pi_num);
                den_K = conv(den_K, pi_den);
                
                details = [details, sprintf('   - PI controller with Kp = %.4f, Ki = %.4f\n', Kp, Ki)];
                
            case 'PD'
                Kp = target_gain;
                
                % Calculate Kd based on bandwidth and plant characteristics
                if isVeryUnstable
                    % For very unstable plants, increase derivative action
                    Kd = Kp / bandwidth * 3;
                    details = [details, '   - Increasing derivative action for highly unstable plant\n'];
                elseif hasCloseRHPZeros
                    % For plants with RHP zeros, reduce derivative action
                    Kd = Kp / bandwidth * 0.5;
                    details = [details, '   - Reducing derivative action due to RHP zeros\n'];
                else
                    % Standard calculation
                    Kd = Kp / bandwidth;
                end
                
                % Create PD controller with filtering
                pd_num = [Kd, Kp];
                pd_den = [epsilon*Kd, 1];
                
                % Combine with numerator and denominator
                num_K = conv(num_K, pd_num);
                den_K = conv(den_K, pd_den);
                
                details = [details, sprintf('   - PD controller with Kp = %.4f, Kd = %.4f, epsilon = %.4f\n', Kp, Kd, epsilon)];
                
            case 'PID'
                Kp = target_gain;
                
                % Calculate Ki based on bandwidth and plant characteristics
                if plantInfo.hasIntegrator
                    % For plants with integrator, use very small Ki
                    Ki = Kp * bandwidth / 50;
                    details = [details, '   - Plant has integrator - using very small Ki\n'];
                elseif isVeryUnstable
                    % For very unstable plants, reduce integral action
                    Ki = Kp * bandwidth / 30;
                    details = [details, '   - Reducing integral action for highly unstable plant\n'];
                elseif hasCloseRHPZeros
                    % For plants with RHP zeros, reduce integral action
                    Ki = Kp * bandwidth / 20;
                    details = [details, '   - Reducing integral action due to RHP zeros\n'];
                else
                    % Standard calculation
                    Ki = Kp * bandwidth / 10;
                end
                
                % Calculate Kd based on bandwidth and plant characteristics
                if isVeryUnstable
                    % For very unstable plants, increase derivative action
                    Kd = Kp / bandwidth * 3;
                    details = [details, '   - Increasing derivative action for highly unstable plant\n'];
                elseif hasCloseRHPZeros
                    % For plants with RHP zeros, reduce derivative action
                    Kd = Kp / bandwidth * 0.5;
                    details = [details, '   - Reducing derivative action due to RHP zeros\n'];
                else
                    % Standard calculation
                    Kd = Kp / bandwidth * 2;
                end
                
                % Create PID controller with filtering
                pid_num = [Kd, Kp, Ki];
                pid_den = [epsilon*Kd, 1, 0];
                
                % Combine with numerator and denominator
                num_K = conv(num_K, pid_num);
                den_K = conv(den_K, pid_den);
                
                details = [details, sprintf('   - PID controller with Kp = %.4f, Ki = %.4f, Kd = %.4f, epsilon = %.4f\n', Kp, Ki, Kd, epsilon)];
                
            otherwise
                error('Unsupported controller structure: %s', structure);
        end
        
        % Apply numerical conditioning for high-order controllers
        [num_K, den_K] = conditionTransferFunction(num_K, den_K);
        
        % Create the final controller
        K = tf(num_K, den_K);
        
        % Step 5: Verify closed-loop stability
        try
            T = feedback(G*K, 1);
            cl_poles = pole(T);
            is_stable = all(real(cl_poles) < 0);
            
            details = [details, '\n7. Stability Check:\n'];
            
            if is_stable
                details = [details, '   - Closed-loop system is stable!\n'];
                
                % Calculate stability margins
                [Gm, Pm, Wcg, Wcp] = margin(G*K);
                
                if ~isempty(Pm) && ~isnan(Pm) && ~isinf(Pm)
                    details = [details, sprintf('   - Phase margin: %.2f degrees at %.3f rad/s\n', Pm, Wcp)];
                end
                
                if ~isempty(Gm) && ~isnan(Gm) && ~isinf(Gm)
                    details = [details, sprintf('   - Gain margin: %.2f dB at %.3f rad/s\n', 20*log10(Gm), Wcg)];
                end
                
                % Calculate and display approximate closed-loop bandwidth
                try
                    [mag, phase, w] = bode(T);
                    mag = squeeze(mag);
                    idx = find(mag < 0.707, 1, 'first');
                    if ~isempty(idx) && idx > 1
                        actual_BW = w(idx-1);
                        details = [details, sprintf('   - Closed-loop bandwidth: %.3f rad/s\n', actual_BW)];
                    end
                catch
                    details = [details, '   - Could not calculate closed-loop bandwidth\n'];
                end
                
                % Calculate and display approximate complementary sensitivity peak
                try
                    [Mt, wMt] = getPeakGain(T);
                    [Ms, wMs] = getPeakGain(feedback(1, G*K));
                    
                    details = [details, sprintf('   - Complementary sensitivity peak: %.2f at %.3f rad/s\n', Mt, wMt)];
                    details = [details, sprintf('   - Sensitivity peak: %.2f at %.3f rad/s\n', Ms, wMs)];
                    
                    % Add robustness assessment
                    if Mt < 1.2 && Ms < 1.5
                        details = [details, '   - Controller has excellent robustness properties\n'];
                    elseif Mt < 1.5 && Ms < 2
                        details = [details, '   - Controller has good robustness properties\n'];
                    elseif Mt < 2 && Ms < 3
                        details = [details, '   - Controller has acceptable robustness properties\n'];
                    else
                        details = [details, '   - Controller has limited robustness properties\n'];
                    end
                catch
                    details = [details, '   - Could not calculate sensitivity peaks\n'];
                end
            else
                % If unstable, try gain adjustments
                details = [details, '   - Closed-loop system is unstable with initial design\n'];
                details = [details, '   - Attempting gain adjustments to stabilize...\n'];
                
                [num, den] = tfdata(K, 'v');
                is_stable = false;
                
                % Try reducing gain until stable
                for scale = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001]
                    K_test = tf(num * scale, den);
                    T_test = feedback(G * K_test, 1);
                    
                    if all(real(pole(T_test)) < 0)
                        K = K_test;
                        is_stable = true;
                        details = [details, sprintf('   - System stabilized with gain factor: %.4f\n', scale)];
                        
                        % Calculate stability margins
                        [Gm, Pm, Wcg, Wcp] = margin(G*K);
                        
                        if ~isempty(Pm) && ~isnan(Pm) && ~isinf(Pm)
                            details = [details, sprintf('   - Phase margin: %.2f degrees at %.3f rad/s\n', Pm, Wcp)];
                        end
                        
                        if ~isempty(Gm) && ~isnan(Gm) && ~isinf(Gm)
                            details = [details, sprintf('   - Gain margin: %.2f dB at %.3f rad/s\n', 20*log10(Gm), Wcg)];
                        end
                        
                        break;
                    end
                end
                
                if ~is_stable
                    details = [details, '   - WARNING: Could not stabilize the system with gain adjustment\n'];
                    details = [details, '   - Trying more aggressive compensation approach...\n'];
                    
                    % Last resort: create a more aggressive compensation
                    [K_last, details_last] = designEmergencyController(G, structure, options, plantInfo);
                    
                    % Test if this controller stabilizes the system
                    T_last = feedback(G*K_last, 1);
                    
                    if all(real(pole(T_last)) < 0)
                        K = K_last;
                        details = [details, '   - Emergency controller design successful!\n'];
                        details = [details, details_last];
                    else
                        details = [details, '   - WARNING: All compensation attempts failed\n'];
                        details = [details, '   - Manual tuning is recommended for this system\n'];
                    end
                end
            end
        catch ME
            details = [details, '\nError during stability check: ' ME.message '\n'];
        end
        
    catch ME
        % Handle design errors
        details = [details, '\nError during compensation controller design: ' ME.message '\n'];
        
        % Create a simple default controller as fallback
        if strcmp(structure, 'P')
            K = tf(0.1, 1);
        elseif strcmp(structure, 'PI')
            K = tf([0.1, 0.01], [1, 0]);
        elseif strcmp(structure, 'PD')
            K = tf([0.1, 0.1], [0.01, 1]);
        else % PID
            K = tf([0.1, 0.1, 0.01], [0.01, 1, 0]);
        end
        
        details = [details, 'Using very conservative controller as fallback.\n'];
    end
    
    % Verify controller stability
    try
        [num, den] = tfdata(K, 'v');
        K_poles = roots(den);
        if any(real(K_poles) > 0)
            details = [details, '\nWARNING: Controller has unstable poles. Applying stabilization.\n'];
            
            % Stabilize controller poles by reflecting unstable poles
            for i = 1:length(K_poles)
                if real(K_poles(i)) > 0
                    K_poles(i) = -real(K_poles(i)) + imag(K_poles(i))*1i;
                end
            end
            
            % Create new controller with stabilized poles
            den_stable = poly(K_poles);
            K = tf(num, den_stable);
            details = [details, 'Controller poles have been stabilized.\n'];
        end
    catch
        % If pole stabilization fails, leave controller as is
    end
    
    % Add final summary
    details = [details, '\nSummary:\n'];
    details = [details, 'The compensation controller design method directly addresses problematic'];
    details = [details, ' plant dynamics by replacing them with more desirable characteristics.'];
    if plantInfo.isUnstable
        details = [details, ' For this unstable plant, poles in the right half-plane were'];
        details = [details, ' replaced with stable ones while preserving similar frequency characteristics.'];
    end
    if plantInfo.hasRHPZeros
        details = [details, ' Right half-plane zeros were addressed with low-pass filtering.'];
    end
    
    % Add controller structure
    [num, den] = tfdata(K, 'v');
    details = [details, sprintf('\n\nFinal Controller K(s):\n')];
    details = [details, sprintf('Numerator: [%s]\n', mat2str(num, 5))];
    details = [details, sprintf('Denominator: [%s]\n', mat2str(den, 5))];
end

% Helper functions

function [isVeryUnstable, isHighlyOscillatory, hasCloseRHPZeros, plantCategory] = analyzePlantDifficulty(plantInfo)
    % Analyze plant to identify particularly challenging characteristics
    
    % Initialize flags
    isVeryUnstable = false;
    isHighlyOscillatory = false;
    hasCloseRHPZeros = false;
    
    % Check for multiple unstable poles or high instability
    p = plantInfo.poles;
    unstable_poles = p(real(p) > 0);
    
    if length(unstable_poles) > 1
        isVeryUnstable = true;  % Multiple unstable poles
    elseif length(unstable_poles) == 1 && real(unstable_poles) > 5
        isVeryUnstable = true;  % Single but strongly unstable pole
    end
    
    % Check for highly oscillatory modes
    oscillatory_poles = p(abs(imag(p)) > abs(real(p)));
    if ~isempty(oscillatory_poles)
        isHighlyOscillatory = true;
    end
    
    % Check for RHP zeros close to imaginary axis
    z = plantInfo.zeros;
    rhp_zeros = z(real(z) > 0);
    
    if ~isempty(rhp_zeros) && min(real(rhp_zeros)) < 0.5
        hasCloseRHPZeros = true;  % RHP zero close to imaginary axis
    end
    
    % Determine plant category for design guidance
    if isVeryUnstable && hasCloseRHPZeros
        plantCategory = 'Very Difficult - Highly unstable with RHP zeros';
    elseif isVeryUnstable
        plantCategory = 'Difficult - Highly unstable';
    elseif hasCloseRHPZeros && plantInfo.isUnstable
        plantCategory = 'Difficult - Unstable with RHP zeros';
    elseif hasCloseRHPZeros
        plantCategory = 'Challenging - RHP zeros';
    elseif plantInfo.isUnstable
        plantCategory = 'Moderately Difficult - Unstable';
    elseif isHighlyOscillatory
        plantCategory = 'Moderately Difficult - Highly oscillatory';
    elseif plantInfo.isHighOrder
        plantCategory = 'Moderately Challenging - High order';
    else
        plantCategory = 'Standard';
    end
end

function poles = findPoorlyDampedPoles(poles_G)
    % Identify poles with damping ratio less than 0.3
    poorly_damped = [];
    
    for i = 1:length(poles_G)
        pole_i = poles_G(i);
        
        if real(pole_i) < 0 && imag(pole_i) ~= 0
            magnitude = abs(pole_i);
            damping = -real(pole_i) / magnitude;
            
            if damping < 0.3
                poorly_damped = [poorly_damped, pole_i];
            end
        end
    end
    
    poles = poorly_damped;
end

function [damping, sort_idx] = getPoleDamping(poles)
    % Calculate damping ratios for complex poles and sort
    damping = zeros(size(poles));
    
    for i = 1:length(poles)
        if imag(poles(i)) ~= 0
            magnitude = abs(poles(i));
            damping(i) = -real(poles(i)) / magnitude;
        else
            damping(i) = 1;  % Real poles have damping ratio 1
        end
    end
    
    [damping, sort_idx] = sort(damping);  % Sort by damping (ascending)
end

function bw_limit = getUnstableBandwidthLimit(plantInfo)
    % Calculate safe bandwidth limit for unstable plants
    p = plantInfo.poles;
    unstable_poles = p(real(p) > 0);
    
    if isempty(unstable_poles)
        bw_limit = Inf;
        return;
    end
    
    % For unstable poles, bandwidth should typically be 2-3x the real part
    % But for multiple unstable poles, we need to be more conservative
    if length(unstable_poles) > 1
        % Most unstable pole determines limit
        max_real = max(real(unstable_poles));
        bw_limit = max_real * (4 - min(3, length(unstable_poles)));  % Reduce factor as poles increase
    else
        % For a single unstable pole
        max_real = real(unstable_poles);
        bw_limit = max_real * 3;  % Standard factor
    end
    
    % Apply upper bound for very unstable systems
    bw_limit = min(bw_limit, 10);
end

function bw_limit = getRHPZeroBandwidthLimit(plantInfo)
    % Calculate bandwidth limitation due to RHP zeros
    z = plantInfo.zeros;
    rhp_zeros = z(real(z) > 0);
    
    if isempty(rhp_zeros)
        bw_limit = Inf;
        return;
    end
    
    % Find zero with smallest real part (closest to imaginary axis)
    min_real = min(real(rhp_zeros));
    
    % Bandwidth should typically be < 0.5x the real part of the RHP zero
    bw_limit = min_real * 0.5;
end

function target_gain = calculateTargetGain(G, plantInfo, bandwidth)
    % Calculate a suitable target gain for achieving desired bandwidth
    
    % Use frequency response at desired bandwidth
    try
        [mag, ~] = bode(G, bandwidth);
        mag = squeeze(mag);
        
        % Target gain to achieve 0dB at bandwidth
        target_gain = 1 / mag;
        
        % Apply safety factors for different plant types
        if plantInfo.isUnstable
            if length(plantInfo.poles(real(plantInfo.poles) > 0)) > 1
                % Multiple unstable poles - be very conservative
                target_gain = target_gain * 0.3;
            else
                % Single unstable pole - be moderately conservative
                target_gain = target_gain * 0.5;
            end
        elseif plantInfo.hasRHPZeros
            % Non-minimum phase plants need some conservatism
            target_gain = target_gain * 0.7;
        elseif plantInfo.isHighOrder
            % High-order plants benefit from gain margin
            target_gain = target_gain * 0.8;
        end
        
        % Bound the gain to reasonable values
        target_gain = min(max(target_gain, 0.01), 100);
        
    catch
        % Fallback if frequency response fails
        if plantInfo.isUnstable
            target_gain = 0.5;  % Conservative gain for unstable plants
        else
            target_gain = 1.0;  % Default gain for stable plants
        end
    end
end

function value = getOption(options, field, default)
    % Safely extract option with default value
    if isfield(options, field)
        value = options.(field);
    else
        value = default;
    end
end

function [num, den] = conditionTransferFunction(num, den)
    % Apply numerical conditioning to transfer function coefficients
    
    % Scale coefficients to avoid numerical issues
    scale_factor = 1;
    
    % Check for large coefficient values
    max_coeff = max(max(abs(num)), max(abs(den)));
    if max_coeff > 1e6
        scale_factor = 1e6 / max_coeff;
        num = num * scale_factor;
    end
    
    % Check for small coefficient values
    min_num = min(abs(num(abs(num) > 0)));
    min_den = min(abs(den(abs(den) > 0)));
    min_coeff = min(min_num, min_den);
    
    if min_coeff < 1e-6
        scale_factor = 1e-6 / min_coeff;
        num = num * scale_factor;
    end
    
    % Remove very small coefficients (numerical noise)
    num(abs(num) < 1e-10 * max(abs(num))) = 0;
    den(abs(den) < 1e-10 * max(abs(den))) = 0;
end

function [peakgain, wpeak] = getPeakGain(sys)
    % Calculate peak gain of frequency response
    w = logspace(-3, 3, 500);
    [mag, ~] = bode(sys, w);
    mag = squeeze(mag);
    
    [peakgain, idx] = max(mag);
    wpeak = w(idx);
end