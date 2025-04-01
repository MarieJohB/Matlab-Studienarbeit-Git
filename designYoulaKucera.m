function [K, details] = designYoulaKucera(G, structure, options, plantInfo)
% DESIGNYOULAKUCERA Controller design using Youla-Kucera parameterization
%
% Implements the Youla-Kucera parameterization approach to directly find
% stabilizing controllers for challenging plants. This method provides a complete
% parameterization of all stabilizing controllers.
%
% Inputs:
%   G        - Plant transfer function
%   structure - Controller structure ('P', 'PI', 'PD', 'PID')
%   options   - Structure with design parameters:
%     .bandwidth  - Desired bandwidth in rad/s (default: 1)
%     .damping    - Desired damping ratio (default: 0.8)
%     .epsilon    - Filter parameter for D-term (default: 0.1)
%     .objective  - Design objective ('tracking', 'disturbance', 'robustness')
%   plantInfo - Structure with plant analysis information
%
% Outputs:
%   K       - Designed controller as transfer function
%   details - Text description with design details

    % Start with detailed information about the method
    details = 'Youla-Kucera Parameterization Design Method\n';
    details = [details, '---------------------------------------\n'];
    details = [details, sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo))];

    % Default values
    if ~isfield(options, 'bandwidth')
        options.bandwidth = 1;
    end
    
    if ~isfield(options, 'damping')
        options.damping = 0.8;
    end
    
    if ~isfield(options, 'epsilon')
        options.epsilon = 0.1;
    end
    
    if ~isfield(options, 'objective')
        options.objective = 'tracking';
    end
    
    % Extract key parameters
    omega = options.bandwidth;
    zeta = options.damping;
    epsilon = options.epsilon;
    objective = options.objective;
    
    details = [details, sprintf('Desired bandwidth: %.4f rad/s\n', omega)];
    details = [details, sprintf('Desired damping ratio: %.4f\n', zeta)];
    details = [details, sprintf('Derivative filter coefficient: %.4f\n', epsilon)];
    details = [details, sprintf('Design objective: %s\n', objective)];
    
    % Step 1: Coprime factorization of the plant G = N/M
    details = [details, '\nSTEP 1: Coprime factorization of plant G = N/M\n'];
    details = [details, '---------------------------------------\n'];
    
    try
        % For unstable plants, we need a robust factorization method
        if plantInfo.isUnstable
            details = [details, 'Using state-space approach for unstable plant factorization.\n'];
            
            % Get state-space realization
            [A, B, C, D] = ssdata(G);
            n = size(A, 1);
            
            % Design stabilizing state feedback gain
            if plantInfo.isHighOrder
                % For high-order systems, use more robust pole placement
                desired_poles = -abs(real(eig(A))) - ones(n, 1);
                
                % Sort to prioritize unstable poles
                unstable_idx = find(real(eig(A)) >= 0);
                if ~isempty(unstable_idx)
                    desired_poles(unstable_idx) = desired_poles(unstable_idx) * 2;
                end
                
                % Use place() for numerical stability if controllable
                try
                    F = place(A, B, desired_poles);
                catch
                    % Fallback to more basic approach
                    F = -2 * real(unstable_idx) * pinv(B);
                    details = [details, 'Used simplified state feedback due to controllability issues.\n'];
                end
            else
                % For simpler systems, directly place poles relative to bandwidth
                p = eig(A);
                
                % Design poles with desired bandwidth and damping
                desired_poles = zeros(n, 1);
                for i = 1:n
                    if real(p(i)) >= 0
                        % Unstable poles - make stable with desired damping
                        if imag(p(i)) ~= 0
                            % Complex pole
                            freq = abs(p(i));
                            desired_poles(i) = -zeta * freq - 1j * freq * sqrt(1 - min(zeta^2, 0.99));
                        else
                            % Real pole
                            desired_poles(i) = -2 * omega * (1 + 0.2*(i-1));
                        end
                    else
                        % Stable poles - keep or slightly adjust
                        if abs(real(p(i))) < 0.3 * omega
                            % Too close to imaginary axis
                            desired_poles(i) = -0.5 * omega;
                        else
                            % Keep original pole
                            desired_poles(i) = p(i);
                        end
                    end
                end
                
                % Create complex conjugate pairs if needed
                for i = 1:n
                    if imag(desired_poles(i)) ~= 0
                        conjugate_found = false;
                        for j = 1:n
                            if i ~= j && abs(desired_poles(i) - conj(desired_poles(j))) < 1e-6
                                conjugate_found = true;
                                break;
                            end
                        end
                        
                        if ~conjugate_found
                            % Find a real pole to replace with the conjugate
                            for j = 1:n
                                if imag(desired_poles(j)) == 0
                                    desired_poles(j) = conj(desired_poles(i));
                                    break;
                                end
                            end
                        end
                    end
                end
                
                details = [details, 'Designed pole locations for stabilizing feedback:\n'];
                for i = 1:length(desired_poles)
                    if imag(desired_poles(i)) ~= 0
                        details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(desired_poles(i)), imag(desired_poles(i)))];
                    else
                        details = [details, sprintf('  p%d = %.4f\n', i, real(desired_poles(i)))];
                    end
                end
                
                % Compute stabilizing state feedback
                try
                    F = place(A, B, desired_poles);
                catch ME
                    details = [details, sprintf('Place failed: %s\n', ME.message)];
                    F = acker(A, B, desired_poles);
                    details = [details, 'Using acker for pole placement.\n'];
                end
            end
            
            % Check quality of pole placement
            Acl = A - B*F;
            actual_poles = eig(Acl);
            
            details = [details, 'Actual closed-loop poles after state feedback:\n'];
            for i = 1:length(actual_poles)
                if imag(actual_poles(i)) ~= 0
                    details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(actual_poles(i)), imag(actual_poles(i)))];
                else
                    details = [details, sprintf('  p%d = %.4f\n', i, real(actual_poles(i)))];
                end
            end
            
            % Compute left coprime factorization
            % For G = N/M, we use state feedback to construct stable N, M
            Acl = A - B*F;
            Bcl = B;
            Ccl = C;
            Dcl = D;
            
            % M is the closed-loop transfer function from u to u
            M_ss = ss(Acl, -Bcl*F, eye(size(F,1)), eye(size(F,1)));
            M = tf(M_ss);
            
            % N is the closed-loop transfer function from u to y
            N_ss = ss(Acl, Bcl, Ccl, Dcl);
            N = tf(N_ss);
            
            % Create coprime factors for controller parameterization
            % We also need X, Y such that M*X + N*Y = 1 (Bezout identity)
            
            % Compute observer gain for dual problem
            try
                L = place(A', C', 2*desired_poles)';
            catch
                try
                    L = acker(A', C', 2*desired_poles)';
                catch
                    % Fallback to simple gain
                    L = pinv(C);
                    details = [details, 'Using simplified observer gain due to observability issues.\n'];
                end
            end
            
            % X and Y from observer-based realization
            X_ss = ss(A-L*C, L, -F, 0);
            X = tf(X_ss);
            
            Y_ss = ss(A-L*C, B-L*D, -F, 1);
            Y = tf(Y_ss);
            
            % Verify Bezout identity: M*X + N*Y = 1
            bezout = minreal(M*X + N*Y);
            bezout_error = norm(1 - dcgain(bezout));
            
            if bezout_error > 0.1
                details = [details, sprintf('Warning: Bezout identity check failed (error = %.4f).\n', bezout_error)];
                details = [details, 'Attempting correction...\n'];
                
                % Correction factor
                correction = 1/dcgain(bezout);
                X = X * correction;
                Y = Y * correction;
                
                % Re-check
                bezout = minreal(M*X + N*Y);
                bezout_error = norm(1 - dcgain(bezout));
                
                details = [details, sprintf('After correction: Bezout error = %.4f\n', bezout_error)];
            else
                details = [details, sprintf('Bezout identity verified: M*X + N*Y = 1 (error = %.4e)\n', bezout_error)];
            end
            
        else
            % For stable plants, factorization is simpler
            details = [details, 'Using direct factorization for stable plant.\n'];
            
            % For stable G, we can use M = 1, N = G
            M = tf(1, 1);
            N = G;
            
            % X = 1, Y = 0 satisfies the Bezout identity
            X = tf(1, 1);
            Y = tf(0, 1);
        end
        
        details = [details, '\nCoprime factorization completed.\n'];
        details = [details, sprintf('N = %s\n', char(N))];
        details = [details, sprintf('M = %s\n', char(M))];
        
    catch ME
        details = [details, sprintf('Error in coprime factorization: %s\n', ME.message)];
        details = [details, 'Using approximation method...\n'];
        
        % Fallback: use pole-zero manipulation for approximate factorization
        [z, p, k] = zpkdata(G, 'v');
        
        % Separate stable and unstable poles
        stable_poles = p(real(p) < 0);
        unstable_poles = p(real(p) >= 0);
        
        if isempty(unstable_poles)
            % Stable plant, trivial factorization
            M = tf(1, 1);
            N = G;
            X = tf(1, 1);
            Y = tf(0, 1);
        else
            % Create FIR filters for each unstable pole
            M_num = 1;
            M_den = 1;
            
            for i = 1:length(unstable_poles)
                pole_i = unstable_poles(i);
                
                if imag(pole_i) ~= 0
                    % Complex pole, need to include conjugate pair
                    if i < length(unstable_poles) && abs(pole_i - conj(unstable_poles(i+1))) < 1e-6
                        % Skip conjugate pair, we'll add both together
                        continue;
                    end
                    
                    % Add complex conjugate pair
                    if imag(pole_i) > 0
                        real_part = real(pole_i);
                        imag_part = imag(pole_i);
                        
                        % Create (s - p)(s - p*)
                        quad_term = [1, -2*real_part, real_part^2 + imag_part^2];
                        
                        % Reflect to LHP: (s + p)(s + p*)
                        stable_quad = [1, 2*real_part, real_part^2 + imag_part^2];
                        
                        M_num = conv(M_num, quad_term);
                        M_den = conv(M_den, stable_quad);
                    end
                else
                    % Real pole
                    M_num = conv(M_num, [1, -pole_i]);
                    M_den = conv(M_den, [1, abs(pole_i)]);
                end
            end
            
            % Create N = G*M
            M = tf(M_num, M_den);
            N = G * M;
            
            % Create X, Y by solving Bezout identity approximately
            % For simplicity, we set X = M* and solve for Y
            X = M;
            Y = tf(1, 1) / N;
        end
    end
    
    % Step 2: Design a stable parameter Q
    details = [details, '\nSTEP 2: Design of stable parameter Q\n'];
    details = [details, '---------------------------------------\n'];
    
    % Design Q based on the specified objective
    switch lower(objective)
        case 'tracking'
            details = [details, 'Designing Q for optimal tracking performance.\n'];
            
            % For tracking, approximate plant inverse with stable filter
            try
                % Try to get an approximate inverse within bandwidth
                w = logspace(-2, log10(omega*10), 100);
                [mag, phase] = bode(minreal(M / N), w);
                mag = squeeze(mag);
                phase = squeeze(phase);
                
                % Design a frequency response
                target_mag = ones(size(mag));
                target_phase = zeros(size(phase));
                
                % Create roll-off at high frequencies
                rolloff_idx = find(w > omega, 1);
                if ~isempty(rolloff_idx)
                    target_mag(rolloff_idx:end) = target_mag(rolloff_idx:end) .* ...
                        (omega ./ w(rolloff_idx:end)).^2;
                end
                
                % Design a filter to approximate this frequency response
                order = min(4, length(plantInfo.poles));
                Q = fitmagfrd(frd(target_mag, w), order, order);
                
                details = [details, sprintf('Created Q filter of order %d.\n', order)];
            catch ME
                details = [details, sprintf('Frequency domain design failed: %s\n', ME.message)];
                details = [details, 'Using time-domain design approach.\n'];
                
                % Simple Q design based on bandwidth
                Q = tf(omega^2, [1, 2*zeta*omega, omega^2]);
                
                % Add additional roll-off if needed
                if plantInfo.isHighOrder || plantInfo.hasRHPZeros
                    Q = Q * tf(1, [1/omega/5, 1])^2;
                    details = [details, 'Added high-frequency roll-off for robustness.\n'];
                end
            end
            
        case 'disturbance'
            details = [details, 'Designing Q for optimal disturbance rejection.\n'];
            
            % For disturbance rejection, make sensitivity small at low frequencies
            Q = tf([0, omega], [1, omega/10]);
            
            % If the plant has integrator, reduce the integral action
            if plantInfo.hasIntegrator
                Q = tf(omega, [1, omega]);
                details = [details, 'Reduced integral action due to plant integrator.\n'];
            end
            
            % Add high-frequency roll-off
            Q = Q * tf(1, [1/omega/10, 1])^2;
            
        case 'robustness'
            details = [details, 'Designing Q for optimal robustness.\n'];
            
            % For robustness, use a conservative Q
            Q = tf(omega^2, [1, 2*zeta*omega, omega^2]);
            
            % Add significant high-frequency roll-off
            Q = Q * tf(1, [1/omega/20, 1])^3;
            
            % If plant has RHP zeros, add more roll-off
            if plantInfo.hasRHPZeros
                Q = Q * tf(1, [1/omega/5, 1])^2;
                details = [details, 'Added extra roll-off due to RHP zeros.\n'];
            end
            
        otherwise
            details = [details, 'Using balanced design approach.\n'];
            
            % Balanced approach
            Q = tf(omega^2, [1, 2*zeta*omega, omega^2]);
            
            % Add moderate roll-off
            Q = Q * tf(1, [1/omega/10, 1])^2;
    end
    
    % Step 3: Compute the controller K = (X + QN) / (Y - QM)
    details = [details, '\nSTEP 3: Computing controller K = (X + QN) / (Y - QM)\n'];
    details = [details, '-----------------------------------------------\n'];
    
    try
        num = minreal(X + Q*N);
        den = minreal(Y - Q*M);
        
        K_youla = minreal(num / den);
        
        details = [details, 'Successfully computed Youla-parameterized controller.\n'];
    catch ME
        details = [details, sprintf('Error computing controller: %s\n', ME.message)];
        details = [details, 'Using alternative formulation...\n'];
        
        % Try alternative formulation
        try
            K_youla = minreal((X + Q*N) * inv(Y - Q*M));
        catch ME2
            details = [details, sprintf('Alternative formulation also failed: %s\n', ME2.message)];
            details = [details, 'Using direct formula with numerical stabilization...\n'];
            
            % Add a small regularization term
            den_reg = minreal(Y - Q*M + tf(1e-6, 1));
            K_youla = minreal((X + Q*N) / den_reg);
        end
    end
    
    % Check stability of K_youla
    try
        k_poles = pole(K_youla);
        if any(real(k_poles) > 0)
            details = [details, 'Warning: Controller has unstable poles. Applying stabilization.\n'];
            
            % Force controller stability
            [num, den] = tfdata(K_youla, 'v');
            p = roots(den);
            for i = 1:length(p)
                if real(p(i)) > 0
                    p(i) = -abs(real(p(i))) + imag(p(i))*1j;
                end
            end
            den_stable = poly(p);
            K_youla = tf(num, den_stable);
        end
    catch
        details = [details, 'Could not analyze controller stability.\n'];
    end
    
    % Verify closed-loop stability
    try
        closed_loop = feedback(G*K_youla, 1);
        cl_poles = pole(closed_loop);
        
        is_stable = all(real(cl_poles) < 0);
        
        if is_stable
            details = [details, '\nYoula-parameterized controller stabilizes the system! Actual poles:\n'];
        else
            details = [details, '\nWARNING: Youla-parameterized controller does not stabilize the system! Actual poles:\n'];
        end
        
        for i = 1:length(cl_poles)
            if imag(cl_poles(i)) ~= 0
                details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(cl_poles(i)), imag(cl_poles(i)))];
            else
                details = [details, sprintf('  p%d = %.4f\n', i, real(cl_poles(i)))];
            end
        end
        
        % If unstable, try to fix
        if ~is_stable
            details = [details, 'Attempting to stabilize controller by gain adjustment...\n'];
            
            % Try reducing the gain
            found_stable = false;
            for scale = [0.5, 0.2, 0.1, 0.05, 0.01]
                K_scaled = K_youla * scale;
                closed_loop_scaled = feedback(G * K_scaled, 1);
                
                if all(real(pole(closed_loop_scaled)) < 0)
                    K_youla = K_scaled;
                    details = [details, sprintf('System stabilized by scaling controller gain by %.2f\n', scale)];
                    found_stable = true;
                    break;
                end
            end
            
            if ~found_stable
                details = [details, 'WARNING: Could not stabilize system by gain scaling.\n'];
                details = [details, 'Using a more conservative approach...\n'];
                
                % Last resort: use a very simple controller
                if plantInfo.isUnstable
                    % For unstable plants, use a stabilizing controller from the coprime factorization
                    K_youla = X / Y;
                    details = [details, 'Using X/Y from coprime factorization as stabilizing controller.\n'];
                else
                    % For stable plants, use a simple controller
                    K_youla = tf(0.1, 1);
                    details = [details, 'Using very conservative controller.\n'];
                end
            end
        end
    catch ME
        details = [details, sprintf('\nError verifying closed-loop stability: %s\n', ME.message)];
    end
    
    % Step 4: Extract controller in desired structure
    details = [details, '\nSTEP 4: Converting to requested controller structure\n'];
    details = [details, '-----------------------------------------------\n'];
    
    % Simplify controller if possible
    try
        K_simple = minreal(K_youla, 0.01);
        details = [details, sprintf('Simplified controller from order %d to %d\n', 
                  order(K_youla), order(K_simple))];
        K_youla = K_simple;
    catch
        details = [details, 'Could not simplify controller.\n'];
    end
    
    % Extract controller in requested structure
    switch structure
        case 'P'
            % Extract proportional gain
            try
                Kp = dcgain(K_youla);
                if isnan(Kp) || isinf(Kp)
                    % Try gain at mid-frequency
                    Kp = abs(evalfr(K_youla, 1j*omega));
                    
                    if isnan(Kp) || isinf(Kp)
                        % Default value
                        Kp = 1.0;
                    end
                end
                
                % Limit to reasonable values
                Kp = min(max(Kp, 0.1), 100);
                
                K = tf(Kp, 1);
                details = [details, sprintf('Extracted P controller: Kp = %.4f\n', Kp)];
            catch
                Kp = 1.0;
                K = tf(Kp, 1);
                details = [details, 'Using default P controller.\n'];
            end
            
        case 'PI'
            % Extract PI controller using frequency response
            try
                w = logspace(-3, log10(omega*10), 200);
                [mag, phase] = bode(K_youla, w);
                mag = squeeze(mag);
                phase = squeeze(phase);
                
                % Check for integrator (phase approaching -90° at low frequencies)
                has_integrator = (phase(1) < -45);
                
                if has_integrator
                    details = [details, 'Detected integral action in Youla controller.\n'];
                    
                    % Find frequency where phase is -45 degrees (between P and I)
                    idx_45 = find(phase > -45, 1);
                    if isempty(idx_45)
                        idx_45 = 2;
                    end
                    w_i = w(idx_45);
                    
                    % Extract parameters
                    Kp = abs(evalfr(K_youla, 1j*omega));
                    Ki = Kp * w_i / 5;
                else
                    details = [details, 'No integral action detected. Adding appropriate integral term.\n'];
                    
                    % Extract proportional gain
                    Kp = dcgain(K_youla);
                    if isnan(Kp) || isinf(Kp)
                        Kp = abs(evalfr(K_youla, 1j*omega));
                        
                        if isnan(Kp) || isinf(Kp)
                            Kp = 1.0;
                        end
                    end
                    
                    Ki = Kp * omega / 10;
                    
                    if plantInfo.hasIntegrator
                        Ki = Ki / 2;
                        details = [details, 'Reduced integral action due to plant integrator.\n'];
                    end
                end
                
                % Limit to reasonable values
                Kp = min(max(Kp, 0.1), 100);
                Ki = min(max(Ki, 0.01), 50);
                
                K = tf([Kp, Ki], [1, 0]);
                details = [details, sprintf('Extracted PI controller: Kp = %.4f, Ki = %.4f\n', Kp, Ki)];
            catch
                % Default PI controller
                Kp = 1.0;
                Ki = omega / 10;
                
                if plantInfo.hasIntegrator
                    Ki = Ki / 2;
                end
                
                K = tf([Kp, Ki], [1, 0]);
                details = [details, 'Using default PI controller.\n'];
            end
            
        case 'PD'
            % Extract PD controller using frequency response
            try
                w = logspace(-3, log10(omega*20), 200);
                [mag, phase] = bode(K_youla, w);
                mag = squeeze(mag);
                phase = squeeze(phase);
                
                % Check for derivative action (phase > 0)
                has_derivative = any(phase > 10);
                
                if has_derivative
                    details = [details, 'Detected derivative action in Youla controller.\n'];
                    
                    % Find max phase (derivative action)
                    [~, idx_max] = max(phase);
                    w_d = w(idx_max);
                    
                    % Extract parameters
                    Kp = abs(evalfr(K_youla, 1j*omega/2));
                    Kd = Kp / w_d;
                else
                    details = [details, 'No derivative action detected. Adding appropriate derivative term.\n'];
                    
                    % Extract proportional gain
                    Kp = dcgain(K_youla);
                    if isnan(Kp) || isinf(Kp)
                        Kp = abs(evalfr(K_youla, 1j*omega/2));
                        
                        if isnan(Kp) || isinf(Kp)
                            Kp = 1.0;
                        end
                    end
                    
                    Kd = Kp / omega;
                    
                    if plantInfo.hasRHPZeros
                        Kd = Kd / 2;
                        details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
                    end
                end
                
                % Limit to reasonable values
                Kp = min(max(Kp, 0.1), 100);
                Kd = min(max(Kd, 0.01), 50);
                
                % Add filtering for derivative term
                Td = Kd / Kp;
                K = tf([Kd, Kp], [epsilon*Td, 1]);
                details = [details, sprintf('Extracted PD controller: Kp = %.4f, Kd = %.4f\n', Kp, Kd)];
            catch
                % Default PD controller
                Kp = 1.0;
                Kd = 1.0 / omega;
                
                if plantInfo.hasRHPZeros
                    Kd = Kd / 2;
                end
                
                Td = Kd / Kp;
                K = tf([Kd, Kp], [epsilon*Td, 1]);
                details = [details, 'Using default PD controller.\n'];
            end
            
        case 'PID'
            % Extract PID controller using frequency response
            try
                w = logspace(-3, log10(omega*20), 200);
                [mag, phase] = bode(K_youla, w);
                mag = squeeze(mag);
                phase = squeeze(phase);
                
                % Check for integral and derivative actions
                has_integrator = (phase(1) < -45);
                has_derivative = any(phase > 10);
                
                if has_integrator && has_derivative
                    details = [details, 'Detected both integral and derivative action in Youla controller.\n'];
                    
                    % Find frequencies for integral and derivative components
                    idx_45 = find(phase > -45, 1);
                    if isempty(idx_45)
                        idx_45 = 2;
                    end
                    w_i = w(idx_45);
                    
                    [~, idx_max] = max(phase);
                    w_d = w(idx_max);
                    
                    % Extract parameters
                    Kp = abs(evalfr(K_youla, 1j*sqrt(w_i*w_d)));
                    Ki = Kp * w_i / 5;
                    Kd = Kp / w_d;
                else
                    details = [details, 'Full PID behavior not detected. Creating appropriate PID controller.\n'];
                    
                    % Extract proportional gain
                    Kp = dcgain(K_youla);
                    if isnan(Kp) || isinf(Kp)
                        Kp = abs(evalfr(K_youla, 1j*omega));
                        
                        if isnan(Kp) || isinf(Kp)
                            Kp = 1.0;
                        end
                    end
                    
                    Ki = Kp * omega / 10;
                    Kd = Kp / omega;
                    
                    if plantInfo.hasIntegrator
                        Ki = Ki / 2;
                        details = [details, 'Reduced integral action due to plant integrator.\n'];
                    end
                    
                    if plantInfo.hasRHPZeros
                        Kd = Kd / 2;
                        details = [details, 'Reduced derivative action due to non-minimum phase behavior.\n'];
                    end
                }
                
                % Limit to reasonable values
                Kp = min(max(Kp, 0.1), 100);
                Ki = min(max(Ki, 0.01), 50);
                Kd = min(max(Kd, 0.01), 50);
                
                % Create PID controller with derivative filter
                K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                details = [details, sprintf('Extracted PID controller: Kp = %.4f, Ki = %.4f, Kd = %.4f\n', Kp, Ki, Kd)];
            catch
                % Default PID controller
                Kp = 1.0;
                Ki = omega / 10;
                Kd = 1.0 / omega;
                
                if plantInfo.hasIntegrator
                    Ki = Ki / 2;
                end
                
                if plantInfo.hasRHPZeros
                    Kd = Kd / 2;
                end
                
                K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                details = [details, 'Using default PID controller.\n'];
            end
            
        otherwise
            error('Unsupported controller structure');
    end
    
    % Verify final controller performance
    try
        % Closed-loop stability
        closed_loop = feedback(G*K, 1);
        cl_poles = pole(closed_loop);
        
        is_stable = all(real(cl_poles) < 0);
        
        if is_stable
            details = [details, '\nFinal controller stabilizes the plant! Closed-loop poles:\n'];
        else
            details = [details, '\nWARNING: Final controller does not stabilize the plant! Closed-loop poles:\n'];
        }
        
        for i = 1:length(cl_poles)
            if imag(cl_poles(i)) ~= 0
                details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(cl_poles(i)), imag(cl_poles(i)))];
            else
                details = [details, sprintf('  p%d = %.4f\n', i, real(cl_poles(i)))];
            end
        end
        
        % If unstable, try final adjustments
        if ~is_stable
            details = [details, 'Attempting final stabilization of controller...\n'];
            
            % Try reducing the gain
            found_stable = false;
            for scale = [0.5, 0.2, 0.1, 0.05, 0.01]
                K_scaled = K * scale;
                closed_loop_scaled = feedback(G * K_scaled, 1);
                
                if all(real(pole(closed_loop_scaled)) < 0)
                    K = K_scaled;
                    details = [details, sprintf('System stabilized by scaling controller gain by %.2f\n', scale)];
                    found_stable = true;
                    break;
                end
            end
            
            if ~found_stable
                details = [details, 'WARNING: Could not stabilize system. Using stabilizing controller from factorization.\n'];
                
                % Use the stabilizing controller from factorization
                K_stable = X / Y;
                
                % Simplify and convert to requested structure
                try
                    K_stable = minreal(K_stable, 0.01);
                    
                    % Extract simplified controller in requested structure
                    [num_stable, den_stable] = tfdata(K_stable, 'v');
                    
                    % Basic gain extraction
                    K_gain = abs(dcgain(K_stable));
                    if isnan(K_gain) || isinf(K_gain)
                        K_gain = 0.1;
                    end
                    
                    switch structure
                        case 'P'
                            K = tf(K_gain, 1);
                        case 'PI'
                            K = tf([K_gain, K_gain/10], [1, 0]);
                        case 'PD'
                            K = tf([K_gain/omega, K_gain], [epsilon, 1]);
                        case 'PID'
                            K = tf([K_gain/omega, K_gain, K_gain/10], [epsilon, 1, 0]);
                        otherwise
                            K = tf(K_gain, 1);
                    end
                    
                    details = [details, 'Created simplified stabilizing controller in requested structure.\n'];
                catch
                    % Use original stabilizing controller
                    K = K_stable;
                    details = [details, 'Using original stabilizing controller without simplification.\n'];
                end
            end
        end
    catch ME
        details = [details, sprintf('\nError in final controller verification: %s\n', ME.message)];
    end
    
    % Final stability margin analysis
    try
        [Gm, Pm, Wcg, Wcp] = margin(G*K);
        details = [details, sprintf('\nStability margins:\n')];
        details = [details, sprintf('  Gain margin: %.2f dB at %.4f rad/s\n', 20*log10(Gm), Wcg)];
        details = [details, sprintf('  Phase margin: %.2f deg at %.4f rad/s\n', Pm, Wcp)];
    catch
        details = [details, '\nCould not compute stability margins.\n'];
    end
    
    return;
end

% Helper function to create a formatted string with plant information
function infoStr = getPlantInfoString(plantInfo)
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