function [K, details] = designIMC(G, structure, options, plantInfo)
    % Enhanced Internal Model Control method with improved handling of unstable systems
    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    
    % Extract needed parameters
    epsilon = options.epsilon;
    settlingTime = options.settlingTime;
    
    details = [details, sprintf('Target Settling Time: %.2f s\n', settlingTime)];
    
    % For unstable plants, first pre-stabilize
    if plantInfo.isUnstable
        details = [details, 'Plant is unstable. Using pre-stabilization before IMC method.\n'];
        
        % Create a robust stabilizing controller for the unstable plant
        K_stab = designRobustStabilizingController(G, plantInfo);
        
        % Create a stabilized version of the plant
        G_stab = feedback(G * K_stab, 1);
        G_for_design = G_stab;
        isPreStabilized = true;
    else
        G_for_design = G;
        isPreStabilized = false;
    }
    
    % Apply the IMC method
    try
        % Separate the plant into its minimum phase (invertible) and non-minimum phase parts
        [z, p, k] = zpkdata(G_for_design, 'v');
        
        % Identify non-minimum phase zeros (RHP zeros)
        rhp_zero_idx = find(real(z) > 0);
        
        if ~isempty(rhp_zero_idx)
            % Non-minimum phase plant
            details = [details, 'Plant has RHP zeros. Applying factorization for non-minimum phase system.\n'];
            
            % Replace RHP zeros with their mirror images to create invertible portion
            z_plus = z;
            z_minus = z;
            
            for i = 1:length(rhp_zero_idx)
                idx = rhp_zero_idx(i);
                z_plus(idx) = -z(idx);  % Mirror across imaginary axis
                z_minus(idx) = [];      % Remove from invertible portion
            }
            
            % Create invertible and non-invertible parts
            G_plus = zpk(z_plus, p, k);
            G_minus = zpk(z_minus, [], 1);
        else
            % Minimum phase plant
            G_plus = G_for_design;
            G_minus = 1;
        }
        
        % Calculate time constant based on settling time
        % For 2% settling time, use 4*tau = settlingTime
        tau = settlingTime / 4;
        
        % Create IMC filter
        % For different controller orders:
        switch structure
            case 'P'
                % First-order filter for P controller
                F = tf(1, [tau, 1]);
                
            case 'PI'
                % Second-order filter for PI controller
                F = tf(1, [tau^2, 2*tau, 1]);
                
            case 'PD'
                % First-order filter with extra zero for PD controller
                F = tf([tau/2, 1], [tau, 1]);
                
            case 'PID'
                % Second-order filter with extra zero for PID controller
                F = tf([tau^2/4, tau/2, 1], [tau^2, 2*tau, 1]);
                
            otherwise
                error('Unsupported controller structure for IMC method');
        end
        
        % Generate IMC controller
        Q = minreal(F / G_plus);
        
        % Convert IMC controller to standard feedback controller
        K_imc = minreal(Q / (1 - G_for_design * Q));
        
        % Extract appropriate controller structure
        switch structure
            case 'P'
                % Extract proportional gain
                Kp = dcgain(K_imc);
                
                % Ensure positive gain
                Kp = abs(Kp);
                
                K_struct = tf(Kp, 1);
                details = [details, sprintf('P controller: Kp = %.4f\n', Kp)];
                
            case 'PI'
                % Convert to PI structure
                [num, den] = tfdata(K_imc, 'v');
                
                % Try to fit PI controller
                if length(num) >= 2 && length(den) >= 1
                    % Extract approximate PI parameters
                    if den(1) == 0
                        % Pure integral term exists, use different extraction
                        Kp = num(1) / den(2);
                        Ki = num(2) / den(2);
                    else
                        % Standard case
                        s = tf('s');
                        [Kp, Ki] = getPIparameters(K_imc);
                    end
                    
                    % Ensure positive gains
                    Kp = abs(Kp);
                    Ki = abs(Ki);
                    
                    K_struct = tf([Kp, Ki], [1, 0]);
                    details = [details, sprintf('PI controller: Kp = %.4f, Ki = %.4f\n', Kp, Ki)];
                else
                    % If PI fit not possible, use direct approximation
                    Kp = abs(dcgain(K_imc));
                    Ki = Kp / tau;  % Approximate integral time based on time constant
                    
                    K_struct = tf([Kp, Ki], [1, 0]);
                    details = [details, sprintf('PI controller (approximated): Kp = %.4f, Ki = %.4f\n', Kp, Ki)];
                end
                
            case 'PD'
                % Convert to PD structure
                try
                    % Try to extract PD parameters
                    s = tf('s');
                    [Kp, Kd] = getPDparameters(K_imc);
                    
                    % Ensure positive gains
                    Kp = abs(Kp);
                    Kd = abs(Kd);
                    
                    if Kd < 0.001 * Kp
                        Kd = 0.1 * Kp * tau;  % Minimum derivative action
                    end
                    
                    Td = Kd / Kp;
                    K_struct = tf([Kd, Kp], [epsilon*Td, 1]);
                    details = [details, sprintf('PD controller: Kp = %.4f, Kd = %.4f\n', Kp, Kd)];
                catch
                    % If extraction fails, use approximation
                    Kp = abs(dcgain(K_imc));
                    Kd = 0.1 * Kp * tau;  % Approximate derivative time
                    
                    Td = Kd / Kp;
                    K_struct = tf([Kd, Kp], [epsilon*Td, 1]);
                    details = [details, sprintf('PD controller (approximated): Kp = %.4f, Kd = %.4f\n', Kp, Kd)];
                end
                
            case 'PID'
                % Convert to PID structure
                try
                    % Try to extract PID parameters
                    s = tf('s');
                    [Kp, Ki, Kd] = getPIDparameters(K_imc);
                    
                    % Ensure positive gains
                    Kp = abs(Kp);
                    Ki = abs(Ki);
                    Kd = abs(Kd);
                    
                    if Ki < 0.001 * Kp
                        Ki = 0.1 * Kp / tau;  % Minimum integral action
                    end
                    
                    if Kd < 0.001 * Kp
                        Kd = 0.1 * Kp * tau;  % Minimum derivative action
                    end
                    
                    Td = Kd / Kp;
                    K_struct = tf([Kd, Kp, Ki], [epsilon*Td, 1, 0]);
                    details = [details, sprintf('PID controller: Kp = %.4f, Ki = %.4f, Kd = %.4f\n', Kp, Ki, Kd)];
                catch
                    % If extraction fails, use approximation
                    Kp = abs(dcgain(K_imc));
                    Ki = 0.1 * Kp / tau;  % Approximate integral time
                    Kd = 0.1 * Kp * tau;  % Approximate derivative time
                    
                    Td = Kd / Kp;
                    K_struct = tf([Kd, Kp, Ki], [epsilon*Td, 1, 0]);
                    details = [details, sprintf('PID controller (approximated): Kp = %.4f, Ki = %.4f, Kd = %.4f\n', Kp, Ki, Kd)];
                end
                
            otherwise
                error('Unsupported controller structure for IMC method');
        end
        
        % For pre-stabilized systems, combine with the stabilizing controller
        if isPreStabilized
            K_combined = series(K_struct, K_stab);
            
            try
                % Simplify the combined controller if possible
                K_combined = minreal(K_combined, 0.01);
                details = [details, 'Successfully simplified the combined controller.\n'];
            catch
                details = [details, 'Could not simplify the combined controller.\n'];
            end
            
            K = K_combined;
            details = [details, 'Combined with pre-stabilizing controller for final result.\n'];
        else
            K = K_struct;
        end
        
    catch ME
        details = [details, sprintf('Error in IMC method: %s\n', ME.message)];
        
        % Fallback to conservative controller
        if isPreStabilized
            % If pre-stabilization worked, just use that controller
            K = K_stab;
            details = [details, 'Using pre-stabilizing controller as fallback.\n'];
        else
            % Create a conservative controller based on structure
            switch structure
                case 'P'
                    K = tf(0.1, 1);
                case 'PI'
                    K = tf([0.1, 0.01], [1, 0]);
                case 'PD'
                    K = tf([0.02, 0.1], [epsilon, 1]);
                case 'PID'
                    K = tf([0.02, 0.1, 0.01], [epsilon, 1, 0]);
                otherwise
                    K = tf(0.1, 1);
            end
            details = [details, 'Using default conservative controller as fallback.\n'];
        end
    end
    
    % Verify controller stability
    try
        K_poles = pole(K);
        if any(real(K_poles) > 0)
            details = [details, 'WARNING: Controller has unstable poles. Applying stabilization.\n'];
            
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