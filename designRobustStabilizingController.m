function [K, details] = designRobustStabilizingController(G, plantInfo)
% DESIGNROBUSTSTABILIZINGCONTROLLER Creates a controller that robustly stabilizes difficult plants
%
% Specialized function for creating controllers that stabilize particularly
% challenging plants, especially those with RHP poles, RHP zeros, or both.
%
% Inputs:
%   G         - Plant transfer function
%   plantInfo - Structure with plant analysis information
%
% Outputs:
%   K        - Stabilizing controller
%   details  - String with details about the design approach

    % Initialize details
    details = 'ROBUST STABILIZING CONTROLLER DESIGN\n';
    details = [details, '----------------------------------\n'];
    
    % If plantInfo is not provided, analyze the plant
    if nargin < 2 || isempty(plantInfo)
        plantInfo = analyzePlant(G);
    end
    
    % Log plant characteristics
    details = [details, sprintf('Plant Type: %s\n', plantInfo.type)];
    details = [details, sprintf('Difficulty Level: %s (%d/100)\n', ...
        plantInfo.difficultyLevel, plantInfo.controlDifficulty)];
    
    % Extract key plant information
    [num_G, den_G] = tfdata(G, 'v');
    p = plantInfo.poles;
    z = plantInfo.zeros;
    hasRHPZeros = plantInfo.hasRHPZeros;
    isUnstable = plantInfo.isUnstable;
    
    % If plant is stable, return a simple controller
    if ~isUnstable
        details = [details, 'Plant is already stable. Returning minimal controller.\n'];
        K = tf(1, 1);
        return;
    end
    
    details = [details, sprintf('Unstable poles: %d\n', plantInfo.unstablePoleCount)];
    if hasRHPZeros
        details = [details, sprintf('RHP zeros: %d\n', plantInfo.rhpZeroCount)];
    end
    
    try
        % Strategy 1: Targeted Pole Compensation (for 1-2 unstable poles)
        if plantInfo.unstablePoleCount <= 2 && ~hasRHPZeros
            details = [details, '\nUsing Targeted Pole Compensation strategy.\n'];
            
            % Extract unstable poles
            unstable_poles = p(real(p) > 0);
            
            % Initialize controller TF components
            num_K = 1;
            den_K = 1;
            
            for i = 1:length(unstable_poles)
                pole_i = unstable_poles(i);
                
                if imag(pole_i) ~= 0
                    % Skip complex conjugate pair (handle both at once)
                    if i < length(unstable_poles) && abs(pole_i - conj(unstable_poles(i+1))) < 1e-6
                        continue;
                    end
                    
                    % For complex poles, create a quadratic factor
                    if imag(pole_i) > 0
                        real_part = real(pole_i);
                        imag_part = imag(pole_i);
                        
                        % Create numerator that cancels unstable complex poles
                        quad_term_num = [1, -2*real_part, real_part^2 + imag_part^2];
                        
                        % Create stable replacement with similar frequency
                        stable_real_part = -2 * abs(real_part);
                        quad_term_den = [1, -2*stable_real_part, stable_real_part^2 + imag_part^2];
                        
                        num_K = conv(num_K, quad_term_num);
                        den_K = conv(den_K, quad_term_den);
                        
                        details = [details, sprintf('Compensating complex unstable pole at %.3f+%.3fi\n', real_part, imag_part)];
                    end
                else
                    % For real poles
                    real_pole = real(pole_i);
                    
                    % Create numerator that cancels unstable pole
                    num_K = conv(num_K, [1, -real_pole]);
                    
                    % Create stable replacement with similar time constant
                    stable_pole = -2 * real_pole;
                    den_K = conv(den_K, [1, -stable_pole]);
                    
                    details = [details, sprintf('Compensating real unstable pole at %.3f\n', real_pole)];
                end
            end
            
            % Add appropriate gain
            K_gain = plantInfo.estimatedStabilizingGain;
            if isnan(K_gain) || K_gain == 0
                K_gain = 1;
            end
            
            num_K = num_K * K_gain;
            
            % Create controller transfer function
            K = tf(num_K, den_K);
            
            % Check if this stabilizes the plant
            T = feedback(G*K, 1);
            if all(real(pole(T)) < 0)
                details = [details, 'Targeted compensation successful! Closed-loop system is stable.\n'];
            else
                details = [details, 'Targeted compensation unsuccessful. Trying alternative approach.\n'];
                % Continue to next strategy
            end
        elseif plantInfo.unstablePoleCount <= 2 && hasRHPZeros
            % Strategy 2: Partial pole compensation with RHP zero handling
            details = [details, '\nUsing Partial Pole Compensation with RHP Zero Handling.\n'];
            
            % Extract unstable poles and RHP zeros
            unstable_poles = p(real(p) > 0);
            rhp_zeros = z(real(z) > 0);
            
            % Initialize controller TF components
            num_K = 1;
            den_K = 1;
            
            % Handle case with both RHP poles and zeros
            details = [details, 'Plant has both unstable poles and RHP zeros.\n'];
            
            % Determine if a simple approach can work
            simpleCaseWorks = false;
            
            % For plants with one unstable pole and one RHP zero
            if length(unstable_poles) == 1 && length(rhp_zeros) == 1 && imag(unstable_poles(1)) == 0 && imag(rhp_zeros(1)) == 0
                % Get the real parts
                unstable_pole = real(unstable_poles(1));
                rhp_zero = real(rhp_zeros(1));
                
                % If the zero is "faster" than the pole, we can use a simple approach
                if rhp_zero > unstable_pole * 1.5
                    simpleCaseWorks = true;
                    
                    details = [details, 'RHP zero is faster than unstable pole. Using simplified approach.\n'];
                    
                    % Create a proportional-like controller with high-frequency roll-off
                    K_gain = plantInfo.estimatedStabilizingGain * 1.5; % Higher gain for robustness
                    
                    if isnan(K_gain) || K_gain == 0
                        K_gain = unstable_pole * 2;
                    end
                    
                    % Create controller: K * (s + a) / (s + b)
                    a = unstable_pole * 0.8; % Slightly less than the unstable pole
                    b = unstable_pole * 10;  % Much faster stable pole
                    
                    K = tf(K_gain * [1, a], [1, b]);
                    
                    % Check if this stabilizes the system
                    T = feedback(G*K, 1);
                    if all(real(pole(T)) < 0)
                        details = [details, 'Simple approach successful! Closed-loop system is stable.\n'];
                    else
                        simpleCaseWorks = false;
                        details = [details, 'Simple approach failed. Using more sophisticated method.\n'];
                    end
                end
            end
            
            if ~simpleCaseWorks
                % More complex case requires a more sophisticated approach
                details = [details, 'Using advanced compensation method for RHP zeros.\n'];
                
                % Determine lowest frequency RHP zero
                [min_rhp_zero, idx] = min(real(rhp_zeros));
                
                % Design a controller with bandwidth below the RHP zero
                target_bandwidth = min_rhp_zero * 0.3;
                
                % Create a lead-lag controller
                % Lead term for phase boost
                lead_zero = target_bandwidth * 0.3;
                lead_pole = lead_zero * 0.1;
                lead_term = tf([1, lead_zero], [1, lead_pole]);
                
                % Lag term for low-frequency gain
                lag_zero = target_bandwidth * 0.01;
                lag_pole = lag_zero * 0.1;
                lag_term = tf([1, lag_zero], [1, lag_pole]);
                
                % Combine terms
                K_base = lead_term * lag_term;
                
                % Add gain
                K_gain = plantInfo.estimatedStabilizingGain * 2;
                if isnan(K_gain) || K_gain == 0
                    K_gain = 10;
                end
                
                K = K_base * K_gain;
                
                % Check if this stabilizes the system
                T = feedback(G*K, 1);
                if all(real(pole(T)) < 0)
                    details = [details, 'Advanced compensation successful! Closed-loop system is stable.\n'];
                else
                    details = [details, 'Advanced compensation unsuccessful. Trying final approach.\n'];
                    % Continue to next strategy
                end
            end
        else
            % Strategy 3: High-order compensation for complex cases
            details = [details, '\nUsing High-Order Compensation for Complex Case.\n'];
            details = [details, sprintf('Plant has %d unstable poles and %d RHP zeros.\n', ...
                       plantInfo.unstablePoleCount, plantInfo.rhpZeroCount)];
            
            % For systems with many unstable poles or both RHP poles and zeros
            % We'll use a more aggressive compensation approach
            
            % Determine overall approach based on plant characteristics
            if plantInfo.controlDifficulty >= 80
                details = [details, 'Plant is extremely difficult to control. Using multi-stage approach.\n'];
                
                % For extremely challenging plants, try a sequential approach
                % Stage 1: Design an inner stabilizing loop with minimal performance expectations
                
                % Identify the most unstable pole
                [max_real_pole, idx] = max(real(p));
                
                if imag(p(idx)) == 0
                    % Real dominant unstable pole
                    inner_gain = max_real_pole * 3;
                    inner_zero = max_real_pole * 0.8;
                    inner_pole = max_real_pole * 10;
                    
                    K_inner = tf(inner_gain * [1, inner_zero], [1, inner_pole]);
                else
                    % Complex dominant unstable pole
                    inner_gain = abs(p(idx)) * 3;
                    K_inner = tf(inner_gain, 1);
                end
                
                % Check if inner controller stabilizes
                G_inner = feedback(G * K_inner, 1);
                inner_poles = pole(G_inner);
                
                if any(real(inner_poles) > 0)
                    details = [details, 'Inner loop failed to stabilize. Trying direct approach.\n'];
                    
                    % Fall back to a more direct approach
                    % Use a very aggressive PID-like controller
                    
                    % Estimate appropriate parameters based on plant characteristics
                    Kp = plantInfo.estimatedStabilizingGain * 5;
                    if isnan(Kp) || Kp == 0
                        Kp = 100;  % Very high gain as a last resort
                    end
                    
                    % Add derivative action for phase lead
                    Kd = Kp / max_real_pole;
                    
                    % Add minimal integral action
                    Ki = Kp * max_real_pole * 0.01;
                    
                    % Create PID with heavy filtering
                    epsilon = 0.01;  % Strong filtering
                    K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                    
                    details = [details, 'Created high-gain filtered PID controller as last resort.\n'];
                else
                    details = [details, 'Inner loop successfully stabilized! Building outer loop.\n'];
                    
                    % Now design outer loop for better performance
                    % Use a simple PI controller for the stabilized plant
                    Kp_outer = 1;
                    Ki_outer = 0.1;
                    K_outer = tf([Kp_outer, Ki_outer], [1, 0]);
                    
                    % Combine inner and outer controllers
                    K_combined = series(K_outer, feedback(K_inner, G));
                    
                    % Try to simplify the combined controller
                    try
                        K = minreal(K_combined);
                    catch
                        K = K_combined;
                    end
                    
                    details = [details, 'Successfully created two-stage controller.\n'];
                end
            else
                % For moderately difficult plants, try a direct approach
                details = [details, 'Using direct compensation approach.\n'];
                
                % Create a controller that attempts to stabilize with reasonable gain
                % and appropriate phase characteristics
                
                % Estimate appropriate gain
                K_gain = plantInfo.estimatedStabilizingGain * 2;
                if isnan(K_gain) || K_gain == 0
                    K_gain = 10;  % Default high gain
                end
                
                % Create controller with phase lead at critical frequencies
                avg_unstable_pole = mean(real(p(real(p) > 0)));
                
                % Lead compensator with zero near the average unstable pole
                lead_zero = avg_unstable_pole * 0.8;
                lead_pole = lead_zero * 0.1;
                
                % Add roll-off for high frequencies
                rolloff_pole = avg_unstable_pole * 10;
                
                % Combine components
                num_K = K_gain * conv([1, lead_zero], [1, lead_zero/5]);
                den_K = conv([1, lead_pole], [1, rolloff_pole]);
                
                % If there are RHP zeros, add low-pass filtering
                if hasRHPZeros
                    min_rhp_zero = min(real(rhp_zeros));
                    lp_cutoff = min_rhp_zero * 0.3;
                    den_K = conv(den_K, [1/lp_cutoff^2, sqrt(2)/lp_cutoff, 1]);
                end
                
                K = tf(num_K, den_K);
                
                details = [details, 'Created direct compensation controller.\n'];
            end
        end
        
        % Final check for stabilization
        try
            T = feedback(G*K, 1);
            cl_poles = pole(T);
            is_stable = all(real(cl_poles) < 0);
            
            details = [details, '\nFinal Stability Check:\n'];
            
            if is_stable
                details = [details, 'SUCCESS: Controller successfully stabilizes the plant!\n'];
                
                % Calculate stability margins if possible
                try
                    [Gm, Pm, ~, ~] = margin(G*K);
                    details = [details, sprintf('Gain Margin: %.2f dB\n', 20*log10(Gm))];
                    details = [details, sprintf('Phase Margin: %.2f degrees\n', Pm)];
                catch
                    details = [details, 'Could not calculate stability margins.\n'];
                end
                
                % Calculate sensitivity peak
                try
                    S = feedback(1, G*K);
                    [Ms, ~] = getPeakGain(S);
                    details = [details, sprintf('Sensitivity Peak: %.2f\n', Ms)];
                catch
                    details = [details, 'Could not calculate sensitivity peak.\n'];
                end
            else
                details = [details, 'WARNING: Controller does not stabilize the plant!\n'];
                
                % Last resort: try a very simple high-gain controller
                details = [details, 'Attempting emergency stabilization...\n'];
                
                % Very simple high-gain controller as a last resort
                K = tf(1000, [1, 0]);
                
                % Check if this extreme approach works
                T = feedback(G*K, 1);
                if all(real(pole(T)) < 0)
                    details = [details, 'Emergency stabilization successful with pure integrator.\n'];
                else
                    details = [details, 'All stabilization attempts failed for this extremely challenging plant.\n'];
                    details = [details, 'Consider reformulating the control problem or plant model.\n'];
                end
            end
        catch ME
            details = [details, sprintf('Error in final stability check: %s\n', ME.message)];
        end
        
    catch ME
        % If something went wrong in the design process
        details = [details, sprintf('Error in controller design: %s\n', ME.message)];
        
        % Create a simple default controller as emergency fallback
        K = tf(100, [1, 0]);
        details = [details, 'Created emergency high-gain controller as fallback.\n'];
    end
end

% Helper function to calculate peak gain of a system
function [peakgain, wpeak, w, mag] = getPeakGain(sys)
    % Get frequency range adjusted to system dynamics
    p = pole(sys);
    z = zero(sys);
    
    % Find appropriate frequency range
    if ~isempty(p) && ~isempty(z)
        % Use pole-zero information to set frequency range
        max_freq = max([10 * max(abs(real(p))) + max(abs(imag(p))), ...
                        10 * max(abs(real(z))) + max(abs(imag(z)))]);
    else
        % Default range if pole/zero info not available
        max_freq = 1000;
    end
    
    w = logspace(-3, log10(max_freq), 500);
    
    % Calculate magnitude response
    [mag, ~] = bode(sys, w);
    mag = squeeze(mag);
    
    % Find peak and corresponding frequency
    [peakgain, idx] = max(mag);
    wpeak = w(idx);
end