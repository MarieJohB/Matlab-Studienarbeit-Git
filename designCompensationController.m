function [K, details] = designCompensationController(G, structure, options, plantInfo)
% DESIGNCOMPENSATIONCONTROLLER Design a controller that directly compensates for problematic plant dynamics
% This method is particularly useful for plants with problematic poles or zeros, including unstable systems.
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
    
    % Extract options
    bandwidth = options.bandwidth;
    damping = options.damping;
    epsilon = options.epsilon;
    
    % Default controller parameters
    Kp = 1;
    Ki = 0;
    Kd = 0;
    
    try
        % Get plant in zero-pole-gain form
        [zeros_G, poles_G, gain_G] = zpkdata(G, 'v');
        
        % Step 1: Identify problematic dynamics to compensate
        unstable_poles = poles_G(real(poles_G) > 0);
        rhp_zeros = zeros_G(real(zeros_G) > 0);
        oscillatory_poles = poles_G(abs(imag(poles_G)) > 0.1*abs(real(poles_G)) & real(poles_G) < 0);
        poorly_damped = find_poorly_damped_poles(poles_G);

        details = [details, '1. Problematic Dynamics Identified:\n'];
        if ~isempty(unstable_poles)
            details = [details, sprintf('   - %d unstable pole(s)\n', length(unstable_poles))];
        end
        if ~isempty(rhp_zeros)
            details = [details, sprintf('   - %d RHP zero(s)\n', length(rhp_zeros))];
        end
        if ~isempty(oscillatory_poles)
            details = [details, sprintf('   - %d oscillatory pole(s)\n', length(oscillatory_poles)/2)];
        end
        if ~isempty(poorly_damped)
            details = [details, sprintf('   - %d poorly damped pole(s)\n', length(poorly_damped)/2)];
        end
        
        % Step 2: Design compensator components
        num_K = 1;
        den_K = 1;
        
        % Handle unstable poles - use pole cancellation with stable replacement
        if ~isempty(unstable_poles)
            details = [details, '\n2. Compensating Unstable Poles:\n'];
            
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
                        new_real_part = -abs(real_part) * 2; % Increase damping
                        new_quad_term = [1, -2*new_real_part, new_real_part^2 + imag_part^2];
                        
                        num_K = conv(num_K, quad_term_num);
                        den_K = conv(den_K, new_quad_term);
                        
                        details = [details, sprintf('   - Complex pole at %.3f+%.3fi replaced with %.3f+%.3fi\n', real_part, imag_part, new_real_part, imag_part)];
                    end
                else
                    % For real poles
                    % Create numerator to cancel the unstable pole
                    num_K = conv(num_K, [1, -pole_i]);
                    
                    % Create a stable replacement (negative real part with increased magnitude)
                    stable_pole = -abs(pole_i) * 2;
                    den_K = conv(den_K, [1, -stable_pole]);
                    
                    details = [details, sprintf('   - Real pole at %.3f replaced with %.3f\n', pole_i, stable_pole)];
                end
            end
        end
        
        % Handle poorly damped poles - add damping if needed
        if ~isempty(poorly_damped) && isempty(unstable_poles)
            details = [details, '\n2. Improving Poorly Damped Poles:\n'];
            
            for i = 1:length(poorly_damped)
                pole_i = poorly_damped(i);
                
                % Skip conjugate pairs (we'll handle both at once)
                if i < length(poorly_damped) && abs(pole_i - conj(poorly_damped(i+1))) < 1e-6
                    continue;
                end
                
                if imag(pole_i) > 0
                    % For complex poles, create quadratic terms
                    real_part = real(pole_i);
                    imag_part = imag(pole_i);
                    magnitude = abs(pole_i);
                    
                    % Calculate current damping ratio
                    current_damping = -real_part / magnitude;
                    
                    if current_damping < 0.4  % Only improve damping if it's too low
                        % Create numerator to cancel the poorly damped poles
                        quad_term_num = [1, -2*real_part, real_part^2 + imag_part^2];
                        
                        % Create better damped replacement with same frequency
                        new_damping = min(0.7, current_damping * 2);  % Double the damping but cap at 0.7
                        new_real_part = -new_damping * magnitude;
                        new_quad_term = [1, -2*new_real_part, magnitude^2];
                        
                        num_K = conv(num_K, quad_term_num);
                        den_K = conv(den_K, new_quad_term);
                        
                        details = [details, sprintf('   - Increased damping for pole at %.3f+%.3fi from %.2f to %.2f\n', real_part, imag_part, current_damping, new_damping)];
                    end
                end
            end
        end
        
        % Handle RHP zeros - cannot be cancelled but can be mitigated
        if ~isempty(rhp_zeros)
            details = [details, '\n3. Addressing RHP Zeros:\n'];
            details = [details, '   - RHP zeros cannot be cancelled directly\n'];
            details = [details, '   - Adding conservative low-pass filtering to reduce sensitivity\n'];
            
            % Add low-pass filtering based on the rightmost RHP zero
            [max_real_zero, idx] = max(real(rhp_zeros));
            filter_freq = max_real_zero * 0.5;  % Set filter below the frequency of the RHP zero
            
            % Add low-pass filter
            den_K = conv(den_K, [1/filter_freq^2, sqrt(2)/filter_freq, 1]);
            
            details = [details, sprintf('   - Added second-order low-pass filter with cutoff at %.3f rad/s\n', filter_freq)];
        end
        
        % Step 3: Based on requested structure, add appropriate terms
        details = [details, '\n4. Applying Requested Controller Structure:\n'];
        
        % Calculate a gain based on desired bandwidth
        % Use the complementary sensitivity crossover as a rough approximation
        if plantInfo.isUnstable
            % More conservative for unstable plants
            nominal_gain = bandwidth * 0.5;
        else
            nominal_gain = bandwidth;
        end
        
        % Adjust based on plant DC gain if available
        if ~isnan(plantInfo.dcGain) && ~isinf(plantInfo.dcGain) && plantInfo.dcGain ~= 0
            target_gain = nominal_gain / abs(plantInfo.dcGain);
        else
            target_gain = nominal_gain;
        end
        
        % Apply structure-specific controller elements
        switch structure
            case 'P'
                Kp = target_gain;
                
                details = [details, sprintf('   - P controller with Kp = %.4f\n', Kp)];
                
            case 'PI'
                Kp = target_gain;
                
                % Calculate Ki based on bandwidth and desired phase margin
                Ki = Kp * bandwidth / 10;  % Conservative ratio for stability
                
                % Add integrator
                den_K = conv(den_K, [1, 0]);
                num_K = conv(num_K, [Kp, Ki]);
                
                details = [details, sprintf('   - PI controller with Kp = %.4f, Ki = %.4f\n', Kp, Ki)];
                
            case 'PD'
                Kp = target_gain;
                
                % Calculate Kd based on bandwidth and desired phase margin
                Kd = Kp / bandwidth;
                
                % Add filtered derivative
                num_K = conv(num_K, [Kd, Kp]);
                den_K = conv(den_K, [epsilon*Kd, 1]);
                
                details = [details, sprintf('   - PD controller with Kp = %.4f, Kd = %.4f, epsilon = %.4f\n', Kp, Kd, epsilon)];
                
            case 'PID'
                Kp = target_gain;
                
                % Calculate Ki and Kd based on bandwidth and damping
                Ki = Kp * bandwidth / 10;
                Kd = Kp / bandwidth * 2;
                
                % Add integrator and filtered derivative
                num_K = conv(num_K, [Kd, Kp, Ki]);
                den_K = conv(den_K, [epsilon*Kd, 1, 0]);
                
                details = [details, sprintf('   - PID controller with Kp = %.4f, Ki = %.4f, Kd = %.4f, epsilon = %.4f\n', Kp, Ki, Kd, epsilon)];
                
            otherwise
                error('Unsupported controller structure: %s', structure);
        end
        
        % Create the final controller
        K = tf(num_K, den_K);
        
        % Step 5: Verify closed-loop stability
        try
            T = feedback(G*K, 1);
            cl_poles = pole(T);
            is_stable = all(real(cl_poles) < 0);
            
            details = [details, '\n5. Stability Check:\n'];
            
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
            else
                % If unstable, try gain adjustments
                details = [details, '   - Closed-loop system is unstable with initial design\n'];
                details = [details, '   - Attempting gain adjustments to stabilize...\n'];
                
                [num, den] = tfdata(K, 'v');
                is_stable = false;
                
                % Try reducing gain until stable
                for scale = [0.5, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001]
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

function poorly_damped = find_poorly_damped_poles(poles)
    % Identify poles with damping ratio less than 0.3
    poorly_damped = [];
    
    for i = 1:length(poles)
        pole_i = poles(i);
        
        if real(pole_i) < 0 && imag(pole_i) ~= 0
            magnitude = abs(pole_i);
            damping = -real(pole_i) / magnitude;
            
            if damping < 0.3
                poorly_damped = [poorly_damped, pole_i];
            end
        end
    end
end