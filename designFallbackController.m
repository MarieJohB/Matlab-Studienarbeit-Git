function [K, details] = designFallbackController(G, structure, epsilon, plantInfo)
% DESIGNFALLBACKCONTROLLER Create a robust fallback controller for challenging plants
% 
% This function creates a robust controller that attempts to stabilize difficult
% plants when other design methods fail. It's specifically designed to handle
% highly unstable systems by combining approaches from different advanced methods.
%
% Inputs:
%   G        - Plant transfer function
%   structure - Controller structure ('P', 'PI', 'PD', 'PID')
%   epsilon   - Filter parameter for D-term (default: 0.1)
%   plantInfo - Structure with plant analysis information
%
% Outputs:
%   K       - Designed controller as transfer function
%   details - Text description with design details

    details = sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo));
    details = [details, 'Applying robust fallback controller design.\n'];
    
    % Start with pre-stabilization for highly unstable plants
    if plantInfo.isUnstable
        details = [details, 'Using multi-stage approach for unstable plant.\n'];
        
        % Step 1: Identify unstable poles for targeted stabilization
        p = plantInfo.poles;
        unstable_poles = p(real(p) > 0);
        
        details = [details, sprintf('Detected %d unstable poles:\n', length(unstable_poles))];
        for i = 1:length(unstable_poles)
            if imag(unstable_poles(i)) ~= 0
                details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(unstable_poles(i)), imag(unstable_poles(i)))];
            else
                details = [details, sprintf('  p%d = %.4f\n', i, real(unstable_poles(i)))];
            end
        end
        
        % Step 2: Create a stabilizing controller with direct pole-zero cancellation
        try
            details = [details, 'Attempting direct pole cancellation approach.\n'];
            
            % For multiple unstable poles, use pole cancellation
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
            
            % Ensure proper controller (numpoles >= numzeros)
            if length(K_num) > length(K_den)
                % Add poles far in the LHP for proper controller
                diff = length(K_num) - length(K_den);
                fastest_unstable = max(real(unstable_poles));
                extra_pole = [1, 10*fastest_unstable];
                
                for i = 1:diff
                    K_den = conv(K_den, extra_pole);
                end
            end
            
            % Add integrator if needed for steady-state tracking
            if ~plantInfo.hasIntegrator && strcmpi(structure, 'PI') || strcmpi(structure, 'PID')
                K_den = conv(K_den, [1, 0]);
                details = [details, 'Added integrator for improved tracking.\n'];
            end
            
            % Initial stabilizing gain
            K_gain = min(1.0, 1.0 / length(unstable_poles));
            
            % Basic stabilizing controller
            K_stab = tf(K_gain * K_num, K_den);
            
            % Test stability of pre-stabilization
            closed_loop = feedback(G * K_stab, 1);
            if all(real(pole(closed_loop)) < 0)
                details = [details, 'Direct pole cancellation successful!\n'];
            else
                % If direct cancellation fails, adjust gain
                details = [details, 'Initial pole cancellation not successful. Adjusting gain...\n'];
                
                % Try different gains to find a stabilizing controller
                for gain_scale = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005]
                    K_test = tf((K_gain * gain_scale) * K_num, K_den);
                    closed_loop = feedback(G * K_test, 1);
                    
                    if all(real(pole(closed_loop)) < 0)
                        K_stab = K_test;
                        details = [details, sprintf('Found stabilizing gain at %.4f * original gain.\n', gain_scale)];
                        break;
                    end
                end
            end
            
            % Check final pre-stabilization
            closed_loop = feedback(G * K_stab, 1);
            prestab_successful = all(real(pole(closed_loop)) < 0);
            
            if prestab_successful
                details = [details, 'Pre-stabilization successful!\n'];
                
                % Use the pre-stabilized system for simpler controller design
                G_stab = closed_loop;
                
                % Now design a performance controller for the stabilized plant
                if strcmpi(structure, 'P')
                    % Simple P controller
                    Kp = min(1.0, 1.0 / abs(dcgain(G_stab)));
                    if isnan(Kp) || isinf(Kp)
                        Kp = 0.5;
                    end
                    K_perf = tf(Kp, 1);
                    
                elseif strcmpi(structure, 'PI')
                    % PI controller
                    Kp = min(1.0, 0.5 / abs(dcgain(G_stab)));
                    if isnan(Kp) || isinf(Kp)
                        Kp = 0.5;
                    end
                    
                    Ti = 10; % Conservative integral time
                    Ki = Kp / Ti;
                    K_perf = tf([Kp, Ki], [1, 0]);
                    
                elseif strcmpi(structure, 'PD')
                    % PD controller
                    Kp = min(1.0, 0.5 / abs(dcgain(G_stab)));
                    if isnan(Kp) || isinf(Kp)
                        Kp = 0.5;
                    end
                    
                    Td = 1.0; % Conservative derivative time
                    K_perf = tf([Kp*Td, Kp], [epsilon*Td, 1]);
                    
                elseif strcmpi(structure, 'PID')
                    % PID controller
                    Kp = min(1.0, 0.5 / abs(dcgain(G_stab)));
                    if isnan(Kp) || isinf(Kp)
                        Kp = 0.5;
                    end
                    
                    Ti = 10; % Conservative integral time
                    Td = 1.0; % Conservative derivative time
                    Ki = Kp / Ti;
                    K_perf = tf([Kp*Td, Kp, Ki], [epsilon*Td, 1, 0]);
                else
                    % Default P controller
                    Kp = min(1.0, 0.5 / abs(dcgain(G_stab)));
                    if isnan(Kp) || isinf(Kp)
                        Kp = 0.5;
                    end
                    K_perf = tf(Kp, 1);
                end
                
                % Combine pre-stabilizing and performance controllers
                K = minreal(series(K_perf, K_stab));
                
                % Final stability check
                closed_loop_final = feedback(G * K, 1);
                if all(real(pole(closed_loop_final)) < 0)
                    details = [details, 'Final controller stabilizes the plant.\n'];
                else
                    details = [details, 'Final controller does not stabilize the plant. Reverting to stabilizing controller only.\n'];
                    K = K_stab;
                end
            else
                details = [details, 'Pre-stabilization failed. Trying alternative approach.\n'];
                % If pre-stabilization fails, proceed to state-space approach
                K = [];
            end
        catch ME
            details = [details, sprintf('Pole cancellation approach failed: %s\n', ME.message)];
            % Reset K to empty to trigger state-space approach
            K = [];
        end
        
        % Try state-space approach if pole cancellation failed
        if isempty(K) || ~prestab_successful
            try
                details = [details, 'Attempting state-space stabilization approach.\n'];
                
                % Convert to state-space for more robust control
                try
                    [A, B, C, D] = ssdata(G);
                catch
                    % If conversion fails, try balanced realization
                    G_bal = balreal(G);
                    [A, B, C, D] = ssdata(G_bal);
                end
                
                n = size(A, 1);
                
                % Design poles with sufficient damping
                desired_poles = zeros(n, 1);
                orig_poles = eig(A);
                
                for i = 1:n
                    if real(orig_poles(i)) >= 0
                        % Move unstable poles to LHP with sufficient damping
                        if imag(orig_poles(i)) ~= 0
                            % Complex pole
                            omega_i = abs(orig_poles(i));
                            desired_poles(i) = -omega_i - 1j * imag(orig_poles(i));
                        else
                            % Real pole
                            desired_poles(i) = -2 * abs(orig_poles(i));
                        end
                    else
                        % Keep stable poles
                        desired_poles(i) = orig_poles(i);
                    end
                end
                
                % Make sure complex poles come in conjugate pairs
                for i = 1:n
                    if imag(desired_poles(i)) ~= 0
                        conj_found = false;
                        for j = i+1:n
                            if abs(desired_poles(i) - conj(desired_poles(j))) < 1e-10
                                conj_found = true;
                                break;
                            end
                        end
                        
                        if ~conj_found
                            % Find a real pole to convert to conjugate
                            for j = 1:n
                                if j ~= i && imag(desired_poles(j)) == 0
                                    desired_poles(j) = conj(desired_poles(i));
                                    break;
                                end
                            end
                        end
                    end
                end
                
                % Use place() for state feedback design
                try
                    F = place(A, B, desired_poles);
                catch
                    % If place() fails, try acker() or a simpler approach
                    try
                        F = acker(A, B, desired_poles);
                    catch
                        % Simple gain calculation as fallback
                        F = -sum(real(unstable_poles)) * pinv(B);
                    end
                end
                
                % Create controller from state feedback
                Acl = A - B*F;
                Bcl = B;
                Ccl = -F;
                Dcl = 0;
                
                K_sf = ss(Acl, Bcl, Ccl, Dcl);
                K_sf_tf = tf(K_sf);
                
                % Check if K_sf_tf stabilizes the plant
                closed_loop = feedback(G * K_sf_tf, 1);
                sf_successful = all(real(pole(closed_loop)) < 0);
                
                if sf_successful
                    details = [details, 'State feedback approach successful!\n'];
                    
                    % Convert to the requested structure
                    switch structure
                        case 'P'
                            Kp = abs(dcgain(K_sf_tf));
                            if isnan(Kp) || isinf(Kp) || Kp < 0.1
                                Kp = 0.5;
                            end
                            K = tf(Kp, 1);
                            
                        case 'PI'
                            [num, den] = tfdata(K_sf_tf, 'v');
                            Kp = 0.5;
                            Ti = 10; % Conservative integral time
                            Ki = Kp / Ti;
                            K = tf([Kp, Ki], [1, 0]);
                            
                        case 'PD'
                            [num, den] = tfdata(K_sf_tf, 'v');
                            Kp = 0.5;
                            Td = 1.0; % Conservative derivative time
                            K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
                            
                        case 'PID'
                            [num, den] = tfdata(K_sf_tf, 'v');
                            Kp = 0.5;
                            Ti = 10; % Conservative integral time
                            Td = 1.0; % Conservative derivative time
                            Ki = Kp / Ti;
                            K = tf([Kp*Td, Kp, Ki], [epsilon*Td, 1, 0]);
                            
                        otherwise
                            K = K_sf_tf;
                    end
                    
                    % Final stability check
                    closed_loop_final = feedback(G * K, 1);
                    if all(real(pole(closed_loop_final)) < 0)
                        details = [details, 'Final structured controller stabilizes the plant.\n'];
                    else
                        details = [details, 'Structured controller does not stabilize the plant. Using state feedback controller directly.\n'];
                        K = K_sf_tf;
                    end
                else
                    details = [details, 'State feedback approach failed. Using final fallback method.\n'];
                    K = [];
                end
            catch ME
                details = [details, sprintf('State-space approach failed: %s\n', ME.message)];
                K = [];
            end
        end
        
        % Final fallback for extremely challenging systems
        if isempty(K) || ~(exist('sf_successful', 'var') && sf_successful) && ~(exist('prestab_successful', 'var') && prestab_successful)
            details = [details, 'Using last-resort fallback controller for highly unstable system.\n'];
            
            % Very specialized controller for highly unstable systems
            switch structure
                case 'P'
                    % Focus on stabilizing the dominant unstable poles
                    K = tf(0.01, 1);
                    
                case 'PI'
                    % Very conservative PI with minimal integral action
                    K = tf([0.01, 0.0001], [1, 0]);
                    
                case 'PD'
                    % PD with strong derivative action for unstable systems
                    K = tf([0.1, 0.01], [epsilon*10, 1]);
                    
                case 'PID'
                    % PID with emphasis on derivative action, minimal integral action
                    K = tf([0.1, 0.01, 0.0001], [epsilon*5, 1, 0]);
                    
                otherwise
                    K = tf(0.01, 1);
            end
            
            % Try different gains to find a stabilizing controller
            stabilized = false;
            for scale = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001]
                K_test = K * scale;
                try
                    closed_loop = feedback(G * K_test, 1);
                    cl_poles = pole(closed_loop);
                    
                    if all(real(cl_poles) < 0)
                        K = K_test;
                        details = [details, sprintf('Found stabilizing gain at scale factor %.5f.\n', scale)];
                        stabilized = true;
                        break;
                    end
                catch
                    % Continue to next scale factor if evaluation fails
                    continue;
                end
            end
            
            if ~stabilized
                details = [details, 'WARNING: Could not find a stabilizing controller. Manual tuning required.\n'];
            end
        end
        
    else
        % For stable plants, use a more conventional approach
        details = [details, 'Using conventional approach for stable plant.\n'];
        
        % Create a simple conservative controller based on plant DC gain
        try
            dcg = abs(dcgain(G));
            if isnan(dcg) || isinf(dcg) || dcg == 0
                dcg = 1;
            end
            
            % Conservative values based on structure
            switch structure
                case 'P'
                    Kp = 0.5 / dcg;
                    K = tf(Kp, 1);
                    
                case 'PI'
                    Kp = 0.4 / dcg;
                    Ti = 10;  % Conservative integral time
                    Ki = Kp / Ti;
                    K = tf([Kp, Ki], [1, 0]);
                    
                case 'PD'
                    Kp = 0.5 / dcg;
                    Td = 0.1;  % Conservative derivative time
                    K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
                    
                case 'PID'
                    Kp = 0.4 / dcg;
                    Ti = 10;  % Conservative integral time
                    Td = 0.1;  % Conservative derivative time
                    Ki = Kp / Ti;
                    K = tf([Kp*Td, Kp, Ki], [epsilon*Td, 1, 0]);
                    
                otherwise
                    Kp = 0.5 / dcg;
                    K = tf(Kp, 1);
            end
            
            % Check stability
            closed_loop = feedback(G * K, 1);
            is_stable = all(real(pole(closed_loop)) < 0);
            
            if is_stable
                details = [details, 'Conservative controller stabilizes the plant.\n'];
            else
                details = [details, 'Conservative controller does not stabilize the plant. Reducing gain.\n'];
                
                % Try reducing gain until stable
                stabilized = false;
                for scale = [0.5, 0.1, 0.05, 0.01]
                    K_scaled = K * scale;
                    closed_loop = feedback(G * K_scaled, 1);
                    
                    if all(real(pole(closed_loop)) < 0)
                        K = K_scaled;
                        details = [details, sprintf('Stabilized with gain scale factor %.2f.\n', scale)];
                        stabilized = true;
                        break;
                    end
                end
                
                if ~stabilized
                    details = [details, 'Could not stabilize with gain reduction. Using very conservative controller.\n'];
                    
                    % Ultra-conservative controller
                    switch structure
                        case 'P'
                            K = tf(0.01, 1);
                        case 'PI'
                            K = tf([0.01, 0.001], [1, 0]);
                        case 'PD'
                            K = tf([0.01, 0.005], [epsilon, 1]);
                        case 'PID'
                            K = tf([0.01, 0.005, 0.0005], [epsilon, 1, 0]);
                        otherwise
                            K = tf(0.01, 1);
                    end
                end
            end
        catch ME
            details = [details, sprintf('Error in controller design: %s\n', ME.message)];
            details = [details, 'Using very conservative default controller.\n'];
            
            % Default conservative controller
            switch structure
                case 'P'
                    K = tf(0.01, 1);
                case 'PI'
                    K = tf([0.01, 0.001], [1, 0]);
                case 'PD'
                    K = tf([0.01, 0.005], [epsilon, 1]);
                case 'PID'
                    K = tf([0.01, 0.005, 0.0005], [epsilon, 1, 0]);
                otherwise
                    K = tf(0.01, 1);
            end
        end
    end
    
    % Final verification and adjustment
    try
        closed_loop = feedback(G * K, 1);
        cl_poles = pole(closed_loop);
        is_stable = all(real(cl_poles) < 0);
        
        details = [details, '\nFinal closed-loop stability check: '];
        if is_stable
            details = [details, 'STABLE\n'];
        else
            details = [details, 'UNSTABLE\n'];
            details = [details, 'WARNING: Final controller does not stabilize the system.\n'];
            details = [details, 'Manual tuning is strongly recommended.\n'];
        end
        
        % Show closed-loop poles
        details = [details, 'Closed-loop poles:\n'];
        for i = 1:length(cl_poles)
            if imag(cl_poles(i)) ~= 0
                details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(cl_poles(i)), imag(cl_poles(i)))];
            else
                details = [details, sprintf('  p%d = %.4f\n', i, real(cl_poles(i)))];
            end
        end
    catch ME
        details = [details, sprintf('\nError in final verification: %s\n', ME.message)];
    end
end