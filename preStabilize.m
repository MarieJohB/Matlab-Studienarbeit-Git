function K = preStabilize(G, plantInfo)
    % PRESTABILIZE Create a stabilizing controller for unstable plants
    % Enhanced version with better handling of highly unstable systems
    
    if ~plantInfo.isUnstable
        K = tf(1, 1); % Unity controller for already stable plants
        return;
    end
    
    % Extract unstable poles for targeted stabilization
    p = plantInfo.poles;
    unstable_poles = p(real(p) > 0);
    
    % Approach 1: Direct pole cancellation for multiple unstable poles
    if length(unstable_poles) > 1 || max(real(unstable_poles)) > 5
        % For highly unstable plants with multiple unstable poles
        % Create a stabilizing controller via pole-zero cancellation
        
        K_num = 1;
        K_den = 1;
        
        % Cancel each unstable pole with a zero
        for i = 1:length(unstable_poles)
            pole_i = unstable_poles(i);
            
            if imag(pole_i) ~= 0
                % Skip conjugate pairs, we'll handle them together
                if i < length(unstable_poles) && abs(pole_i - conj(unstable_poles(i+1))) < 1e-6
                    continue;
                end
                
                if imag(pole_i) > 0
                    real_part = real(pole_i);
                    imag_part = imag(pole_i);
                    
                    % Add zeros at unstable pole locations
                    K_num = conv(K_num, [1, -2*real_part, real_part^2 + imag_part^2]);
                    
                    % Add stable poles further in the LHP
                    K_den = conv(K_den, [1, 5*abs(real_part), 10*(real_part^2 + imag_part^2)]);
                end
            else
                % Real pole
                K_num = conv(K_num, [1, -pole_i]);
                K_den = conv(K_den, [1, 5*abs(pole_i)]);
            end
        end
        
        % Add proper dynamics if necessary
        if length(K_num) > length(K_den)
            % Add poles far in the LHP to ensure proper controller
            diff = length(K_num) - length(K_den);
            fastest_unstable = max(real(unstable_poles));
            extra_pole = [1, 10*fastest_unstable];
            
            for i = 1:diff
                K_den = conv(K_den, extra_pole);
            end
        end
        
        % Scale gain conservatively based on degree of instability
        K_gain = 0.1 / length(unstable_poles);
        if max(real(unstable_poles)) > 10
            K_gain = K_gain * 0.1; % Even more conservative for extremely unstable systems
        end
        
        K = tf(K_gain * K_num, K_den);
        
        % Test stability and adjust gain if needed
        try
            closed_loop = feedback(G*K, 1);
            cl_poles = pole(closed_loop);
            
            if any(real(cl_poles) > 0)
                % Try reducing gain until stable
                for scale = [0.5, 0.1, 0.05, 0.01, 0.005, 0.001]
                    K_test = tf(K_gain * scale * K_num, K_den);
                    closed_loop = feedback(G*K_test, 1);
                    
                    if all(real(pole(closed_loop)) < 0)
                        K = K_test;
                        break;
                    end
                end
            end
        catch
            % If feedback analysis fails, continue with original K
        end
        
    else
        % Approach 2: For single unstable pole, use simpler approach
        pole_i = unstable_poles(1);
        
        if imag(pole_i) ~= 0
            % Complex pole - use second-order approach
            real_part = real(pole_i);
            imag_part = imag(pole_i);
            
            % Quadratic zero term to cancel the complex pole pair
            K_num = [1, -2*real_part, real_part^2 + imag_part^2];
            
            % Stable poles further in LHP
            K_den = [1, 4*abs(real_part), 4*(real_part^2 + imag_part^2)];
            
            K = tf(K_num, K_den);
        else
            % Real pole - design a lead controller
            a = real(pole_i);
            K = tf([1, -a], [1, 2*a]);
            
            % Add high-frequency dynamics for better robustness
            if a > 5 % For very unstable poles
                K = K * tf(1, [1/(10*a), 1]);
            end
            
            % Adjust gain if very unstable
            if a > 10
                K = K * 0.1;
            end
        end
        
        % Test stability and adjust if necessary
        try
            closed_loop = feedback(G*K, 1);
            cl_poles = pole(closed_loop);
            
            if any(real(cl_poles) > 0)
                % Try reducing gain
                for scale = [0.5, 0.2, 0.1, 0.05, 0.01]
                    K_test = K * scale;
                    closed_loop = feedback(G*K_test, 1);
                    
                    if all(real(pole(closed_loop)) < 0)
                        K = K_test;
                        break;
                    end
                end
            end
        catch
            % If feedback analysis fails, continue with original K
        end
    end
    
    % Add small derivative action if high-frequency dynamics present
    if plantInfo.isHighOrder
        K = K * tf([0.05, 1], [0.005, 1]);
    end
    
    % Final verification and extreme fallback if needed
    try
        closed_loop = feedback(G*K, 1);
        if any(real(pole(closed_loop)) > 0)
            % Last resort stabilizer - pure gain with very small value
            disp('Warning: Standard pre-stabilization failed. Using ultra-conservative approach.');
            K = tf(0.0001, 1);
        end
    catch
        % If all else fails, use ultra-conservative gain
        K = tf(0.0001, 1);
    end
end