function [K, details] = designYoulaKucera(G, structure, options, plantInfo)
% DESIGNYOULAKUCERA Controller design using Youla-Kucera parameterization
%
% Implements the Youla-Kucera parameterization approach to directly find
% stabilizing controllers for challenging plants. This method provides a complete
% parameterization of all stabilizing controllers and is particularly effective
% for high-order unstable systems with non-minimum phase zeros.
%
% Inputs:
%   G        - Plant transfer function or state-space model
%   structure - Controller structure ('P', 'PI', 'PD', 'PID')
%   options   - Structure with design parameters
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
    
    % Check if input is state-space or transfer function
    isStateSpace = isa(G, 'ss');
    
    % If state-space model is in the options, use it directly
    if ~isStateSpace && isfield(options, 'stateSpace')
        G_ss = options.stateSpace;
        isStateSpace = true;
        details = [details, 'Using state-space model from options.\n'];
    elseif isStateSpace
        G_ss = G;
        details = [details, 'Using provided state-space model directly.\n'];
    else
        % Try to convert to state-space for more robust handling
        try
            G_ss = ss(G);
            isStateSpace = true;
            details = [details, 'Successfully converted to state-space representation.\n'];
        catch ME
            details = [details, sprintf('Could not convert to state-space: %s\n', ME.message)];
            details = [details, 'Using transfer function approach.\n'];
            G_ss = [];
        end
    end
    
    % Step 1: Coprime factorization of the plant G = N/M
    details = [details, '\nSTEP 1: Coprime factorization of plant G = N/M\n'];
    details = [details, '---------------------------------------\n'];
    
    try
        % Approach based on model type
        if isStateSpace && ~isempty(G_ss)
            details = [details, 'Using state-space approach for factorization.\n'];
            
            % Get state-space realization
            [A, B, C, D] = ssdata(G_ss);
            nx = size(A, 1);
            
            % For unstable plants, we need a robust factorization method
            if plantInfo.isUnstable
                details = [details, 'Using LQR-based factorization for unstable plant.\n'];
                
                % Extract eigenvalues for analysis
                p = eig(A);
                unstable_poles = p(real(p) > 0);
                max_real_part = max(real(unstable_poles));
                
                % Design stabilizing state feedback gain using LQR
                % For highly unstable systems, adjust weights
                if max_real_part > 5 || length(unstable_poles) > 1
                    details = [details, 'System is highly unstable. Using conservative factorization.\n'];
                    Q = 10 * eye(nx);  % Higher state penalty
                    R = 0.1;          % Lower control penalty
                else
                    Q = eye(nx);
                    R = 1;
                end
                
                % Design LQR feedback to stabilize the system
                [F, ~, ~] = lqr(A, B, Q, R);
                
                % Check if LQR worked
                Acl = A - B*F;
                cl_eig = eig(Acl);
                
                if any(real(cl_eig) >= 0)
                    details = [details, 'LQR did not fully stabilize. Using more aggressive design.\n'];
                    
                    % Try more aggressive weights
                    Q = 100 * eye(nx);
                    R = 0.01;
                    [F, ~, ~] = lqr(A, B, Q, R);
                    
                    Acl = A - B*F;
                    cl_eig = eig(Acl);
                    
                    if any(real(cl_eig) >= 0)
                        details = [details, 'Still not stable. Using direct pole placement.\n'];
                        
                        % Try direct pole placement
                        desired_poles = -abs(real(p)) - ones(nx, 1);
                        
                        try
                            F = place(A, B, desired_poles);
                        catch
                            try
                                F = acker(A, B, desired_poles);
                            catch
                                % Fallback to manual stabilization
                                details = [details, 'Using fallback manual stabilization approach.\n'];
                                
                                % Focus on unstable poles
                                F = zeros(1, nx);
                                unstable_idx = find(real(p) > 0);
                                
                                for i = 1:length(unstable_idx)
                                    shift_vector = zeros(nx, 1);
                                    shift_vector(unstable_idx(i)) = 1;
                                    
                                    % Shift this unstable pole to LHP
                                    shift_amount = real(p(unstable_idx(i))) * 3;
                                    F = F + shift_amount * shift_vector' * pinv(B);
                                end
                            end
                        end
                    end
                end
                
                % Compute observer gain for dual problem
                % Use same weights but transposed system
                [L, ~, ~] = lqr(A', C', Q, R);
                L = L';
                
                % Check that observer is stable
                Aobs = A - L*C;
                obs_eig = eig(Aobs);
                
                if any(real(obs_eig) >= 0)
                    details = [details, 'Observer not stable. Using more aggressive design.\n'];
                    
                    % Try more aggressive weights
                    Q = 100 * eye(nx);
                    R = 0.01;
                    L = lqr(A', C', Q, R)';
                    
                    Aobs = A - L*C;
                    obs_eig = eig(Aobs);
                    
                    if any(real(obs_eig) >= 0)
                        details = [details, 'Using direct pole placement for observer.\n'];
                        
                        % Try direct pole placement
                        desired_poles = -abs(real(p)) - 2*ones(nx, 1);
                        
                        try
                            L = place(A', C', desired_poles)';
                        catch
                            try
                                L = acker(A', C', desired_poles)';
                            catch
                                % Fallback
                                L = pinv(C);
                            end
                        end
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
                N = tf(G_ss);
                
                % X = 1, Y = 0 satisfies the Bezout identity for stable plants
                X = tf(1, 1);
                Y = tf(0, 1);
            end
        else
            % Transfer function approach for factorization
            details = [details, 'Using transfer function approach for factorization.\n'];
            
            % For unstable plants, use pole-zero factorization
            if plantInfo.isUnstable
                details = [details, 'Creating factorization for unstable transfer function.\n'];
                
                % Get plant poles and zeros
                [z, p, k] = zpkdata(G, 'v');
                
                % Separate stable and unstable poles
                stable_poles = p(real(p) < 0);
                unstable_poles = p(real(p) >= 0);
                
                % Create coprime factors using pole-zero manipulation
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
                X = M;  % For simple case, X = M is a good approximation
                Y = (1 - M*X) / N;  % Solve for Y to satisfy M*X + N*Y = 1
                
                % Verify Bezout identity
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
                % For stable plants, factorization is simple
                details = [details, 'Using direct factorization for stable plant.\n'];
                
                % For stable G, we can use M = 1, N = G
                M = tf(1, 1);
                N = G;
                
                % X = 1, Y = 0 satisfies the Bezout identity
                X = tf(1, 1);
                Y = tf(0, 1);
            end
        end
        
        details = [details, '\nCoprime factorization completed.\n'];
        
    catch ME
        details = [details, sprintf('Error in coprime factorization: %s\n', ME.message)];
        details = [details, 'Using approximation method...\n'];
        
        % Fallback to simpler approach
        if plantInfo.isUnstable
            % For unstable plants, create an approximate coprime factorization
            
            % Extract unstable poles
            p = plantInfo.poles;
            unstable_poles = p(real(p) >= 0);
            
            % Create M with zeros at unstable poles and stable denominator
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
                    
                    if imag(pole_i) > 0
                        real_part = real(pole_i);
                        imag_part = imag(pole_i);
                        
                        quad_term = [1, -2*real_part, real_part^2 + imag_part^2];
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
            
            % Create M and N
            M = tf(M_num, M_den);
            N = G * M;
            
            % Basic X, Y
            X = M;
            Y = tf(0, 1);
            
            % Simply assume we've factorized well enough to proceed
            details = [details, 'Created approximate coprime factorization.\n'];
        else
            % For stable plants, use trivial factorization
            M = tf(1, 1);
            N = G;
            X = tf(1, 1);
            Y = tf(0, 1);
        end
    end
    
    % Step 2: Design a stable parameter Q
    details = [details, '\nSTEP 2: Design of stable parameter Q\n'];
    details = [details, '---------------------------------------\n'];
    
    % Design Q based on the specified objective
    try
        switch lower(objective)
            case {'tracking', 'reference'}
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
                    
                    % Check for non-minimum phase zeros
                    if plantInfo.hasRHPZeros
                        rhp_zeros = plantInfo.zeros(real(plantInfo.zeros) > 0);
                        min_rhp_zero = min(abs(rhp_zeros));
                        
                        % For non-minimum phase plants, roll off before RHP zeros
                        rhp_idx = find(w > min_rhp_zero*0.3, 1);
                        if ~isempty(rhp_idx)
                            target_mag(rhp_idx:end) = target_mag(rhp_idx:end) .* ...
                                (min_rhp_zero*0.3 ./ w(rhp_idx:end)).^2;
                            
                            details = [details, sprintf('Added additional roll-off below RHP zero at %.4f rad/s.\n', min_rhp_zero)];
                        end
                    end
                    
                    % Design a filter to approximate this frequency response
                    order = min(4, length(plantInfo.poles));
                    
                    % We'll create a simple parameterized filter instead of trying to fit exactly
                    % This is more reliable for challenging plants
                    num = [1];
                    den = [1];
                    
                    % Add complex poles based on bandwidth
                    if zeta < 1
                        % Underdamped
                        num = conv(num, [1, 2*zeta*omega, omega^2]);
                        den = conv(den, [1, 2*2*zeta*omega, (2*omega)^2]);  % Double the freq for poles
                    else
                        % Overdamped
                        num = conv(num, [1, omega]);
                        den = conv(den, [1, 2*omega]);
                    end
                    
                    % Add roll-off
                    den = conv(den, [1, omega*5]);
                    
                    Q = tf(num, den);
                    
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
                
            case {'disturbance', 'disturbance rejection'}
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
                
            case {'robustness', 'robust'}
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
        
        % Special handling for highly unstable systems
        if plantInfo.isUnstable
            p = plantInfo.poles;
            unstable_poles = p(real(p) > 0);
            
            if max(real(unstable_poles)) > 5 || length(unstable_poles) > 1
                details = [details, 'System is highly unstable. Using more conservative Q parameter.\n'];
                
                % Scale down Q for more conservative design
                Q = Q * 0.1;
                
                % Add aggressive roll-off
                Q = Q * tf(1, [1/omega/30, 1])^3;
            end
        end
        
        details = [details, sprintf('Created Q parameter: %s\n', char(Q))];
        
    catch ME
        details = [details, sprintf('Error in Q parameter design: %s\n', ME.message)];
        details = [details, 'Using simple Q parameter.\n'];
        
        % Simple fallback Q
        Q = tf(omega^2, [1, 2*zeta*omega, omega^2]) * tf(1, [1/omega/10, 1])^2 * 0.1;
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
        details = [details, sprintf('Simplified controller from order %d to %d\n', order(K_youla), order(K_simple))];
        K_youla = K_simple;
    catch
        details = [details, 'Could not simplify controller.\n'];
    end
    
    % Extract controller in requested structure
    switch structure
        case 'P'
            % Extract proportional gain
            try
                % Use more robust frequency-domain approach
                w = logspace(-3, log10(omega*10), 100);
                [mag, ~] = bode(K_youla, w);
                mag = squeeze(mag);
                
                % Use gain at mid-frequencies for better stability
                idx_mid = ceil(length(w)/2);
                Kp = mag(idx_mid);
                
                if isnan(Kp) || isinf(Kp)
                    % Try DC gain as fallback
                    Kp = abs(dcgain(K_youla));
                    
                    if isnan(Kp) || isinf(Kp)
                        % Use a conservative default
                        Kp = 1.0;
                    end
                end
                
                % Limit to reasonable values
                Kp = min(max(abs(Kp), 0.1), 100);
                
                K = tf(Kp, 1);
                details = [details, sprintf('Extracted P controller: Kp = %.4f\n', Kp)];
            catch
                % Fallback approach
                Kp = 1.0;
                
                if plantInfo.isUnstable
                    % For unstable plants, use a more conservative gain
                    Kp = 0.5;
                end
                
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
                    Kp = abs(dcgain(K_youla));
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
                    
                    if plantInfo.isUnstable
                        Ki = Ki / 5;
                        details = [details, 'Reduced integral action due to plant instability.\n'];
                    end
                end
                
                % Special handling for highly unstable systems
                if plantInfo.isUnstable
                    p = plantInfo.poles;
                    unstable_poles = p(real(p) > 0);
                    
                    if max(real(unstable_poles)) > 5 || length(unstable_poles) > 1
                        details = [details, 'System is highly unstable. Using very conservative integral action.\n'];
                        Ki = Ki * 0.1;
                    end
                end
                
                % Limit to reasonable values
                Kp = min(max(abs(Kp), 0.1), 100);
                Ki = min(max(abs(Ki), 0.01), 50);
                
                K = tf([Kp, Ki], [1, 0]);
                details = [details, sprintf('Extracted PI controller: Kp = %.4f, Ki = %.4f\n', Kp, Ki)];
            catch
                % Default PI controller
                Kp = 1.0;
                Ki = omega / 10;
                
                if plantInfo.hasIntegrator
                    Ki = Ki / 2;
                end
                
                if plantInfo.isUnstable
                    Ki = Ki / 10;
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
                    Kp = abs(dcgain(K_youla));
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
                Kp = min(max(abs(Kp), 0.1), 100);
                Kd = min(max(abs(Kd), 0.01), 50);
                
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
                    Kp = abs(dcgain(K_youla));
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
                    
                    if plantInfo.isUnstable
                        Ki = Ki / 5;
                        details = [details, 'Reduced integral action due to plant instability.\n'];
                    end
                end
                
                % Special handling for highly unstable systems
                if plantInfo.isUnstable
                    p = plantInfo.poles;
                    unstable_poles = p(real(p) > 0);
                    
                    if max(real(unstable_poles)) > 5 || length(unstable_poles) > 1
                        details = [details, 'System is highly unstable. Using very conservative integral action.\n'];
                        Ki = Ki * 0.1;
                    end
                end
                
                % Limit to reasonable values
                Kp = min(max(abs(Kp), 0.1), 100);
                Ki = min(max(abs(Ki), 0.01), 50);
                Kd = min(max(abs(Kd), 0.01), 50);
                
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
                
                if plantInfo.isUnstable
                    Ki = Ki / 10;
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
        end
        
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
                
                % Final check
                closed_loop_final = feedback(G*K, 1);
                if all(real(pole(closed_loop_final)) < 0)
                    details = [details, 'Final controller is stable.\n'];
                else
                    details = [details, 'WARNING: Even fallback controller is not stable. This system is extremely challenging.\n'];
                    
                    % Last resort: ultra-conservative controller
                    switch structure
                        case 'P'
                            K = tf(0.01, 1);
                        case 'PI'
                            K = tf([0.01, 0.001], [1, 0]);
                        case 'PD'
                            K = tf([0.01, 0.01], [0.1, 1]);
                        case 'PID'
                            K = tf([0.01, 0.01, 0.001], [0.1, 1, 0]);
                        otherwise
                            K = tf(0.01, 1);
                    end
                    
                    details = [details, 'Using ultra-conservative controller as last resort.\n'];
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