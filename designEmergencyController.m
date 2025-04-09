function [K, details] = designEmergencyController(G, structure, options, plantInfo)
    % Last resort emergency controller for very challenging systems
    
    details = 'EMERGENCY CONTROLLER DESIGN\n';
    details = [details, '----------------------------\n'];
    
    % Extract plant poles and zeros
    p = plantInfo.poles;
    
    % First attempt: extreme pole shifting
    try
        % Get poles of the plant
        unstable_poles = p(real(p) > 0);
        
        if isempty(unstable_poles)
            unstable_poles = p;  % If no unstable poles, work with all poles
        end
        
        % Create a controller that shifts all poles far into the left half-plane
        num_K = 1;
        den_K = 1;
        
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
                    
                    % Create numerator to cancel the pole
                    quad_term_num = [1, -2*real_part, real_part^2 + imag_part^2];
                    
                    % Create stable replacement far in the LHP
                    stable_real_part = -10 * abs(real_part) - 5;  % Very conservative
                    new_quad_term = [1, -2*stable_real_part, stable_real_part^2 + imag_part^2];
                    
                    num_K = conv(num_K, quad_term_num);
                    den_K = conv(den_K, new_quad_term);
                end
            else
                % For real poles
                num_K = conv(num_K, [1, -pole_i]);
                
                % Create stable replacement far in the LHP
                stable_pole = -10 * abs(pole_i) - 5;  % Very conservative
                den_K = conv(den_K, [1, -stable_pole]);
            end
        end
        
        % Apply very conservative gain
        gain_factor = 0.001;
        num_K = num_K * gain_factor;
        
        % Add filtering based on controller type
        switch structure
            case 'P'
                % Nothing to add
            case 'PI'
                % Add integral action with very low gain
                ki_factor = 0.0001;
                num_K = conv(num_K, [1, ki_factor]);
                den_K = conv(den_K, [1, 0]);
            case 'PD'
                % Add filtered derivative with very low gain
                kd_factor = 0.0001;
                num_K = conv(num_K, [kd_factor, 1]);
                den_K = conv(den_K, [0.1*kd_factor, 1]);
            case 'PID'
                % Add PID structure with very low gains
                ki_factor = 0.0001;
                kd_factor = 0.0001;
                num_K = conv(num_K, [kd_factor, 1, ki_factor]);
                den_K = conv(den_K, [0.1*kd_factor, 1, 0]);
        end
        
        % Create the final controller
        K = tf(num_K, den_K);
        
        details = [details, 'Created very conservative emergency controller with extreme pole shifting\n'];
        
    catch ME
        % If the complex approach fails, create a minimal stabilizing controller
        details = [details, 'Emergency pole shifting failed: ' ME.message '\n'];
        details = [details, 'Creating minimal stabilizing controller\n'];
        
        % Create a minimal controller based on structure
        switch structure
            case 'P'
                K = tf(0.001, 1);
            case 'PI'
                K = tf([0.001, 0.0001], [1, 0]);
            case 'PD'
                K = tf([0.001, 0.001], [0.001, 1]);
            case 'PID'
                K = tf([0.001, 0.001, 0.0001], [0.001, 1, 0]);
            otherwise
                K = tf(0.001, 1);
        end
    end
    
    % Verify stability and apply additional gain reduction if needed
    try
        T = feedback(G*K, 1);
        if any(real(pole(T)) > 0)
            details = [details, 'Initial emergency controller is still unstable\n'];
            details = [details, 'Applying additional gain reduction\n'];
            
            [num, den] = tfdata(K, 'v');
            
            for scale = [0.1, 0.01, 0.001, 0.0001]
                K_test = tf(num * scale, den);
                T_test = feedback(G * K_test, 1);
                
                if all(real(pole(T_test)) < 0)
                    K = K_test;
                    details = [details, sprintf('System stabilized with additional gain reduction factor: %.6f\n', scale)];
                    break;
                end
            end
        else
            details = [details, 'Emergency controller successfully stabilizes the system\n'];
        end
    catch
        details = [details, 'Could not verify stability of emergency controller\n'];
    end
    
    details = [details, 'WARNING: This is an emergency controller designed for stability only!\n'];
    details = [details, 'Performance will be very conservative and likely slow.\n'];
    details = [details, 'Manual refinement is strongly recommended.\n'];
end