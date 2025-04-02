function K_stab = designRobustStabilizingController(G, plantInfo)
    % Design a robust stabilizing controller for unstable plants
    
    % Extract poles for focused stabilization
    p = plantInfo.poles;
    unstable_poles = p(real(p) > 0);
    
    % Different approaches based on number of unstable poles
    if length(unstable_poles) > 2 || max(real(unstable_poles)) > 10
        % For highly unstable or many poles, use careful pole-zero cancellation
        K_num = 1;
        K_den = 1;
        
        % Create a controller that places zeros at unstable pole locations
        for i = 1:length(unstable_poles)
            pole_i = unstable_poles(i);
            
            if imag(pole_i) ~= 0
                % Skip conjugate pairs, we'll add both together
                if i < length(unstable_poles) && abs(pole_i - conj(unstable_poles(i+1))) < 1e-6
                    continue;
                end
                
                if imag(pole_i) > 0
                    real_part = real(pole_i);
                    imag_part = imag(pole_i);
                    
                    % Create zeros at unstable poles
                    quad_term = [1, -2*real_part, real_part^2 + imag_part^2];
                    
                    % Create stable poles at reflected locations with extra damping
                    stable_quad = [1, 4*real_part, 5*(real_part^2 + imag_part^2)];
                    
                    K_num = conv(K_num, quad_term);
                    K_den = conv(K_den, stable_quad);
                end
            else
                % For real poles
                K_num = conv(K_num, [1, -pole_i]);
                K_den = conv(K_den, [1, 3*pole_i]);  // More damping
            end
        end
        
        % Conservative gain for highly unstable systems
        gain_factor = 0.01;
        
        K_stab = tf(gain_factor * K_num, K_den);
    else
        % For simple unstable systems (1-2 poles), use direct approach
        if length(unstable_poles) == 1 && imag(unstable_poles(1)) == 0
            % Single real unstable pole
            p_unstable = real(unstable_poles(1));
            
            % Simple lead controller
            K_stab = tf([1, -p_unstable], [1, p_unstable*3]);
        else
            % Complex or multiple poles
            K_num = 1;
            K_den = 1;
            
            for i = 1:length(unstable_poles)
                pole_i = unstable_poles(i);
                
                if imag(pole_i) ~= 0
                    % Skip conjugate processing
                    if i < length(unstable_poles) && abs(pole_i - conj(unstable_poles(i+1))) < 1e-6
                        continue;
                    end
                    
                    if imag(pole_i) > 0
                        real_part = real(pole_i);
                        imag_part = imag(pole_i);
                        
                        % Quadratic terms for complex conjugate pair
                        K_num = conv(K_num, [1, -2*real_part, real_part^2 + imag_part^2]);
                        K_den = conv(K_den, [1, 2*real_part, 2*(real_part^2 + imag_part^2)]);
                    end
                else
                    % Real pole
                    K_num = conv(K_num, [1, -pole_i]);
                    K_den = conv(K_den, [1, 2*pole_i]);
                end
            end
            
            K_stab = tf(K_num, K_den);
        end
    end
    
    % Verify and adjust the stabilizing controller
    try
        closed_loop = feedback(G*K_stab, 1);
        cl_poles = pole(closed_loop);
        
        if any(real(cl_poles) > 0)
            % If not stable, try scaling the gain down
            for scale = [0.5, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001]
                [num, den] = tfdata(K_stab, 'v');
                K_test = tf(num * scale, den);
                closed_loop = feedback(G*K_test, 1);
                
                if all(real(pole(closed_loop)) < 0)
                    K_stab = K_test;
                    break;
                end
            end
        end
    catch
        % If analysis fails, use conservative gain
        [num, den] = tfdata(K_stab, 'v');
        K_stab = tf(num * 0.001, den);
    end
    
    return;
end