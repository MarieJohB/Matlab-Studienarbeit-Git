function [K, details] = designEnhancedStateFeedback(G, structure, options, plantInfo)
% DESIGNENHANCEDSTATEFEEDBACK Advanced controller design using state feedback with observer
%
% Implements robust state feedback with observer specifically designed for
% unstable and challenging systems. The method combines pole placement with
% state estimation and includes special handling for various challenging plant types.
%
% Inputs:
%   G        - Plant transfer function or state-space model
%   structure - Controller structure ('P', 'PI', 'PD', 'PID')
%   options   - Structure with design parameters:
%     .bandwidth  - Desired bandwidth in rad/s (default: 1)
%     .damping    - Desired damping ratio (default: 0.8)
%     .epsilon    - Filter parameter for D-term (default: 0.1)
%     .robustness - Robustness level: 'Low', 'Medium', 'High' (default: 'Medium')
%   plantInfo - Structure with plant analysis information
%
% Outputs:
%   K       - Designed controller as transfer function
%   details - Text description with design details

    % Start with detailed information about the method
    details = 'Enhanced State Feedback with Observer Design Method\n';
    details = [details, '--------------------------------------------\n'];
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
    
    if ~isfield(options, 'robustness')
        options.robustness = 'Medium';
    end
    
    % Extract key parameters
    omega = options.bandwidth;
    zeta = options.damping;
    epsilon = options.epsilon;
    robustness = options.robustness;
    
    details = [details, sprintf('Desired bandwidth: %.4f rad/s\n', omega)];
    details = [details, sprintf('Desired damping ratio: %.4f\n', zeta)];
    details = [details, sprintf('Derivative filter coefficient: %.4f\n', epsilon)];
    details = [details, sprintf('Robustness level: %s\n', robustness)];
    
    % Detect if G is already in state-space form
    isStateSpace = isa(G, 'ss');
    
    if isStateSpace
        details = [details, 'Using state-space model directly for controller design.\n'];
        sys_ss = G;
    else
        % Convert to state-space format with careful handling of problematic cases
        try
            sys_ss = ss(G);
            details = [details, 'Successfully converted to state-space representation.\n'];
        catch ME
            details = [details, sprintf('Direct state-space conversion failed: %s\n', ME.message)];
            details = [details, 'Attempting alternative conversion approach...\n'];
            
            try
                % Use alternative balancing approach for better numerical conditioning
                sys_ss = balreal(ss(G));
                details = [details, 'Used balanced realization for improved numerical properties.\n'];
            catch ME2
                details = [details, sprintf('Alternative conversion failed: %s\n', ME2.message)];
                details = [details, 'Attempting model order reduction...\n'];
                
                try
                    % Try model reduction if the system is too complex
                    Gred = balred(G, min(8, length(pole(G))));
                    sys_ss = ss(Gred);
                    details = [details, 'Used model reduction to obtain workable state-space model.\n'];
                catch ME3
                    details = [details, 'All conversion attempts failed. Creating approximate model.\n'];
                    
                    % Create approximate model from poles and zeros
                    z = zero(G);
                    p = pole(G);
                    k = 1;
                    
                    try
                        k = dcgain(G);
                        if isinf(k) || isnan(k)
                            k = 1;
                        end
                    catch
                        % Use default gain
                    end
                    
                    % For problematic cases, create a simplified balanced model
                    A = diag(p);
                    B = ones(length(p), 1);
                    C = ones(1, length(p));
                    D = 0;
                    
                    if ~isempty(z)
                        % Try to incorporate zeros
                        C = [ones(1, length(p)-length(z)), zeros(1, length(z))];
                    end
                    
                    sys_ss = ss(A, B, C, D);
                    details = [details, 'Created simplified model from poles and zeros.\n'];
                end
            end
        end
    end
    
    % Extract state-space matrices
    A = sys_ss.A;
    B = sys_ss.B;
    C = sys_ss.C;
    D = sys_ss.D;
    
    % Get system dimensions
    n = size(A, 1); % Number of states
    
    % Check controllability and observability with robust numerical methods
    Co = ctrb(A, B);
    Ob = obsv(A, C);
    
    % Use SVD for more reliable rank calculation
    sv_ctrb = svd(Co);
    sv_obsv = svd(Ob);
    
    % Set threshold for rank determination
    rank_tol = max(size(A)) * eps(norm(A));
    
    controllability_rank = sum(sv_ctrb > rank_tol);
    observability_rank = sum(sv_obsv > rank_tol);
    
    details = [details, sprintf('Controllability check: rank %d out of %d\n', controllability_rank, n)];
    details = [details, sprintf('Observability check: rank %d out of %d\n', observability_rank, n)];
    
    % Handle non-controllable or non-observable systems
    if controllability_rank < n || observability_rank < n
        details = [details, 'System is not fully controllable or observable.\n'];
        details = [details, 'Using modal decomposition to separate controllable/observable subsystems.\n'];
        
        % Add regularization for numerical stability
        A_reg = A + 1e-10 * eye(n);
        [V, D] = eig(A_reg);
        
        % Sort eigenvalues by controllability/observability measure
        ctrb_measure = zeros(n, 1);
        obsv_measure = zeros(n, 1);
        
        for i = 1:n
            ctrb_measure(i) = norm(B' * V(:,i));
            obsv_measure(i) = norm(C * V(:,i));
        end
        
        % Prioritize unstable modes
        eigvals = diag(D);
        is_unstable = real(eigvals) >= 0;
        
        % Create composite measure, prioritizing unstable modes
        total_measure = ctrb_measure .* obsv_measure;
        total_measure(is_unstable) = total_measure(is_unstable) * 10;
        
        [~, idx] = sort(total_measure, 'descend');
        
        % Keep only controllable and observable modes up to the minimum rank
        n_keep = min(controllability_rank, observability_rank);
        idx_keep = idx(1:n_keep);
        
        % Create reduced model
        V_red = V(:, idx_keep);
        A_red = V_red \ A * V_red;
        B_red = V_red \ B;
        C_red = C * V_red;
        
        % Replace original matrices with reduced controllable/observable system
        A = A_red;
        B = B_red;
        C = C_red;
        n = size(A, 1);
        
        details = [details, sprintf('Reduced to controllable and observable subsystem of order %d.\n', n)];
    end
    
    % Design state feedback gain K_sf and observer gain L
    % Adjust pole locations based on stability properties and robustness
    
    % Calculate desired closed-loop poles for state feedback
    % Adjust based on system type and desired performance
    if plantInfo.isUnstable
        details = [details, 'Unstable system detected. Using conservative pole placement.\n'];
        
        % For unstable systems, determine the unstable poles
        p = eig(A);
        unstable_idx = find(real(p) > 0);
        num_unstable = length(unstable_idx);
        
        details = [details, sprintf('Number of unstable poles: %d\n', num_unstable)];
        
        % For unstable systems, place poles more conservatively
        % Move unstable poles to LHP and keep stable poles relatively unchanged
        sf_poles = p;
        
        % Scale bandwidth based on the most unstable pole
        if ~isempty(unstable_idx)
            most_unstable = max(real(p));
            conservative_omega = max(omega, 2 * most_unstable);
            details = [details, sprintf('Adjusted bandwidth to %.4f rad/s based on unstable poles.\n', conservative_omega)];
        else
            conservative_omega = omega;
        end
        
        % Place unstable poles with sufficient damping
        for i = 1:length(unstable_idx)
            idx = unstable_idx(i);
            if imag(p(idx)) ~= 0
                % Complex unstable pole - maintain frequency but with desired damping
                freq = abs(p(idx));
                sf_poles(idx) = -zeta * freq + 1j * freq * sqrt(1 - zeta^2);
                sf_poles(find(conj(p) == p(idx))) = conj(sf_poles(idx));
            else
                % Real unstable pole - move to LHP with desired frequency
                sf_poles(idx) = -conservative_omega * (1 + 0.5 * (i-1));
            end
        end
        
        % Adjust stable poles slightly if needed
        stable_idx = find(real(p) <= 0);
        for i = 1:length(stable_idx)
            idx = stable_idx(i);
            if real(p(idx)) > -0.2 * conservative_omega
                % If pole is too close to imaginary axis, move it further left
                if imag(p(idx)) ~= 0
                    freq = abs(p(idx));
                    sf_poles(idx) = -zeta * freq + 1j * freq * sqrt(1 - zeta^2);
                    sf_poles(find(conj(p) == p(idx))) = conj(sf_poles(idx));
                else
                    sf_poles(idx) = -0.3 * conservative_omega * (1 + 0.2 * (i-1));
                end
            end
        end
    else
        % For stable systems, use standard second-order pole placement
        if zeta < 1
            % Underdamped - complex conjugate poles
            real_part = -zeta * omega;
            imag_part = omega * sqrt(1 - zeta^2);
            dominant_poles = [complex(real_part, imag_part), complex(real_part, -imag_part)];
        else
            % Critically damped or overdamped - real poles
            dominant_poles = [-omega, -omega * (2*zeta - 1)];
        end
        
        % Additional poles further left for higher-order systems
        if n > 2
            extra_poles = -omega * (3:n+1);
            sf_poles = [dominant_poles, extra_poles(1:(n-2))];
        else
            sf_poles = dominant_poles;
        end
    end
    
    details = [details, 'Desired state feedback poles:\n'];
    for i = 1:length(sf_poles)
        if imag(sf_poles(i)) ~= 0
            details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(sf_poles(i)), imag(sf_poles(i)))];
        else
            details = [details, sprintf('  p%d = %.4f\n', i, real(sf_poles(i)))];
        end
    end
    
    % Calculate observer poles - typically faster than closed-loop poles
    obs_factor = 2.5; % Observer poles 2.5x faster than control poles
    
    % Apply robustness setting to observer poles
    switch robustness
        case 'Low'
            obs_factor = 2.0;  % Faster observer, less filtering
        case 'High'
            obs_factor = 4.0;  % Slower observer, more filtering
        otherwise
            obs_factor = 2.5;  % Medium setting
    end
    
    obs_poles = sf_poles * obs_factor;
    
    details = [details, '\nObserver poles (placed at a factor of ' sprintf('%.1f', obs_factor) ' times control poles):\n'];
    for i = 1:length(obs_poles)
        if imag(obs_poles(i)) ~= 0
            details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(obs_poles(i)), imag(obs_poles(i)))];
        else
            details = [details, sprintf('  p%d = %.4f\n', i, real(obs_poles(i)))];
        end
    end
    
    % Design state feedback and observer gains with robust algorithms
    try
        % Try place() function first for better numerical properties
        K_sf = place(A, B, sf_poles);
        details = [details, 'Successfully computed state feedback gain K using place().\n'];
    catch ME
        details = [details, sprintf('Standard place() failed: %s\n', ME.message)];
        
        try
            % Try acker() as an alternative
            K_sf = acker(A, B, sf_poles);
            details = [details, 'Computed state feedback gain K using acker().\n'];
        catch ME2
            details = [details, sprintf('Acker method failed: %s\n', ME2.message)];
            
            % Use a robust manual pole placement as last resort
            details = [details, 'Using robust manual pole placement.\n'];
            
            % Manual pole placement algorithm
            K_sf = zeros(1, n);
            for i = 1:n
                K_sf = K_sf + abs(real(sf_poles(i))) * pinv(B) * (eye(n) - pinv(A - sf_poles(i)*eye(n)) * A);
            end
        end
    end
    
    details = [details, sprintf('State feedback gain K_sf = [%s]\n', mat2str(K_sf, 4))];
    
    % Design observer gain
    try
        % Observer gain is designed using the dual system
        L = place(A', C', obs_poles)';
        details = [details, 'Successfully computed observer gain L using place().\n'];
    catch ME
        details = [details, sprintf('Observer place() failed: %s\n', ME.message)];
        
        try
            % Try acker() as an alternative
            L = acker(A', C', obs_poles)';
            details = [details, 'Computed observer gain L using acker().\n'];
        catch ME2
            details = [details, sprintf('Observer acker method failed: %s\n', ME2.message)];
            
            % Use a robust manual pole placement as last resort
            details = [details, 'Using robust manual observer gain calculation.\n'];
            
            % Calculate observer gain
            L = zeros(n, 1);
            for i = 1:n
                L = L + abs(real(obs_poles(i))) * (eye(n) - pinv(A' - obs_poles(i)*eye(n)) * A') * pinv(C');
            end
        end
    end
    
    details = [details, sprintf('Observer gain L = [%s]\n', mat2str(L, 4))];
    
    % Create full-order observer-based controller
    Ac = A - B*K_sf - L*C;
    Bc = L;
    Cc = -K_sf;
    Dc = 0;
    
    % Convert to transfer function
    controller_ss = ss(Ac, Bc, Cc, Dc);
    controller_tf = tf(controller_ss);
    
    details = [details, '\nCreated full-order observer-based controller.\n'];
    
    % Convert to the requested controller structure
    details = [details, sprintf('\nConverting observer-based controller to %s structure...\n', structure)];
    
    % Simplify controller transfer function for better numerical properties
    try
        simplified_tf = minreal(controller_tf, 0.01);
        details = [details, sprintf('Reduced controller order from %d to %d through minimization.\n',order(controller_tf), order(simplified_tf))];
        controller_tf = simplified_tf;
    catch
        details = [details, 'Could not simplify controller transfer function.\n'];
    end
    
    % Create the simplified controller in the requested structure
    switch structure
        case 'P'
            % Extract proportional term
            Kp = abs(dcgain(controller_tf));
            if isnan(Kp) || isinf(Kp)
                % For integrating controllers, use gain at low frequency
                w_low = omega/100;
                Kp = abs(evalfr(controller_tf, 1j*w_low)) * w_low;
                if isnan(Kp) || isinf(Kp)
                    Kp = abs(K_sf(1)); % Fallback to first element of state feedback gain
                end
            end
            
            % Limit to reasonable values
            Kp = min(max(Kp, 0.1), 100);
            
            K = tf(Kp, 1);
            details = [details, sprintf('Created P controller: Kp = %.4f\n', Kp)];
            
        case 'PI'
            % Fit a PI controller to the low-frequency response
            w = logspace(-3, log10(omega*10), 50);
            [mag, phase] = bode(controller_tf, w);
            mag = squeeze(mag);
            phase = squeeze(phase);
            
            % Check if controller has integral action
            has_integral = (phase(1) < -45);
            
            if has_integral
                % Extract PI parameters that match low-frequency behavior
                Kp = abs(evalfr(controller_tf, 1j*w(end/4)));
                
                % Find frequency where phase is -45 degrees (between P and I)
                idx = find(phase > -45, 1);
                if isempty(idx)
                    idx = 2;
                end
                w_i = w(idx);
                
                Ki = Kp * w_i / 5;
            else
                % If no integral action, create one
                Kp = abs(dcgain(controller_tf));
                if isnan(Kp) || isinf(Kp)
                    Kp = abs(K_sf(1));
                end
                Ki = Kp * omega / 5;
            end
            
            % Limit to reasonable values for stability
            Kp = min(max(Kp, 0.1), 100);
            Ki = min(max(Ki, 0.01), 50);
            
            K = tf([Kp, Ki], [1, 0]);
            details = [details, sprintf('Created PI controller: Kp = %.4f, Ki = %.4f\n', Kp, Ki)];
            
        case 'PD'
            % For PD controller, we need to extract the derivative action
            w = logspace(-3, log10(omega*10), 50);
            [mag, phase] = bode(controller_tf, w);
            mag = squeeze(mag);
            phase = squeeze(phase);
            
            % Check if controller has derivative action
            has_derivative = any(phase > 45);
            
            if has_derivative
                % Find where phase is maximum (derivative action)
                [~, idx] = max(phase);
                w_d = w(idx);
                
                % Extract PD parameters
                Kp = abs(evalfr(controller_tf, 1j*w(end/4)));
                Kd = Kp / w_d;
            else
                % If no derivative action, estimate from state feedback
                Kp = abs(dcgain(controller_tf));
                if isnan(Kp) || isinf(Kp)
                    Kp = abs(K_sf(1));
                end
                Kd = Kp / omega;
            end
            
            % Limit to reasonable values
            Kp = min(max(Kp, 0.1), 100);
            Kd = min(max(Kd, 0.01), 50);
            Td = Kd / Kp;
            
            % Add filtering for derivative term
            K = tf([Kd, Kp], [epsilon*Td, 1]);
            details = [details, sprintf('Created PD controller: Kp = %.4f, Kd = %.4f\n', Kp, Kd)];
            
        case 'PID'
            % For PID, extract all three actions
            w = logspace(-3, log10(omega*10), 100);
            [mag, phase] = bode(controller_tf, w);
            mag = squeeze(mag);
            phase = squeeze(phase);
            
            % Check controller characteristics
            has_integral = (phase(1) < -45);
            has_derivative = any(phase > 45);
            
            if has_integral && has_derivative
                % Full PID controller
                details = [details, 'Controller shows both integral and derivative characteristics.\n'];
                
                % Find derivative frequency (phase peak)
                [~, idx_d] = max(phase);
                w_d = w(idx_d);
                
                % Find integral frequency (where phase crosses -45°)
                idx_i = find(phase > -45, 1);
                if isempty(idx_i)
                    idx_i = 2;
                end
                w_i = w(idx_i);
                
                % Extract parameters using frequency response
                Kp = abs(evalfr(controller_tf, 1j*sqrt(w_i*w_d)));  % Geometric mean frequency
                Ki = Kp * w_i / 5;
                Kd = Kp / w_d;
            else
                % Partial controller - estimate missing components
                Kp = abs(dcgain(controller_tf));
                if isnan(Kp) || isinf(Kp)
                    Kp = abs(K_sf(1));
                end
                
                if has_integral
                    % PI controller with added derivative
                    idx_i = find(phase > -45, 1);
                    if isempty(idx_i)
                        idx_i = 2;
                    end
                    w_i = w(idx_i);
                    Ki = Kp * w_i / 5;
                    Kd = Kp / omega;
                elseif has_derivative
                    % PD controller with added integral
                    [~, idx_d] = max(phase);
                    w_d = w(idx_d);
                    Kd = Kp / w_d;
                    Ki = Kp * omega / 10;
                else
                    % P controller with added I and D
                    Ki = Kp * omega / 10;
                    Kd = Kp / omega;
                end
            end
            
            % Adjust based on plant characteristics
            if plantInfo.hasIntegrator
                Ki = Ki * 0.5;  % Reduce integral action for plants with integrator
            end
            
            if plantInfo.hasRHPZeros
                Kd = Kd * 0.5;  % Reduce derivative action for non-minimum phase plants
            end
            
            % Limit to reasonable values
            Kp = min(max(Kp, 0.1), 100);
            Ki = min(max(Ki, 0.01), 50);
            Kd = min(max(Kd, 0.01), 50);
            
            % Add filtering for derivative term
            K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
            details = [details, sprintf('Created PID controller: Kp = %.4f, Ki = %.4f, Kd = %.4f\n', Kp, Ki, Kd)];
            
        otherwise
            error('Unsupported controller structure');
    end
    
    % Verify closed-loop stability
    try
        if isStateSpace
            G_tf = tf(G);
            closed_loop = feedback(G_tf*K, 1);
        else
            closed_loop = feedback(G*K, 1);
        end
        
        cl_poles = pole(closed_loop);
        
        is_stable = all(real(cl_poles) < 0);
        
        if is_stable
            details = [details, '\nClosed-loop system is stable! Actual poles:\n'];
        else
            details = [details, '\nWARNING: Closed-loop system is UNSTABLE! Actual poles:\n'];
        end
        
        for i = 1:length(cl_poles)
            if imag(cl_poles(i)) ~= 0
                details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(cl_poles(i)), imag(cl_poles(i)))];
            else
                details = [details, sprintf('  p%d = %.4f\n', i, real(cl_poles(i)))];
            end
        end
        
        % If unstable, attempt stabilization
        if ~is_stable
            details = [details, '\nAttempting to stabilize controller...\n'];
            
            % Try different gain scaling
            scales = [0.5, 0.25, 0.1, 0.05, 0.01];
            for scale = scales
                [num, den] = tfdata(K, 'v');
                K_test = tf(num * scale, den);
                
                if isStateSpace
                    closed_loop_test = feedback(G_tf*K_test, 1);
                else
                    closed_loop_test = feedback(G*K_test, 1);
                end
                
                cl_poles_test = pole(closed_loop_test);
                
                if all(real(cl_poles_test) < 0)
                    K = K_test;
                    details = [details, sprintf('Stabilized by scaling gain by %.3f\n', scale)];
                    
                    % Show final stable poles
                    details = [details, 'Final stable closed-loop poles:\n'];
                    for i = 1:length(cl_poles_test)
                        if imag(cl_poles_test(i)) ~= 0
                            details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(cl_poles_test(i)), imag(cl_poles_test(i)))];
                        else
                            details = [details, sprintf('  p%d = %.4f\n', i, real(cl_poles_test(i)))];
                        end
                    end
                    break;
                end
            end
            
            % If still unstable, try structure-specific adjustments
            if isStateSpace
                closed_loop_check = feedback(G_tf*K, 1);
            else
                closed_loop_check = feedback(G*K, 1);
            end
            
            if ~all(real(pole(closed_loop_check)) < 0)
                details = [details, 'Gain scaling alone could not stabilize. Trying structure-specific modifications.\n'];
                
                switch structure
                    case 'P'
                        % For P controller, just reduce gain more aggressively
                        Kp = dcgain(K) * 0.001;
                        K = tf(Kp, 1);
                        
                    case 'PI'
                        % For PI, reduce integral action drastically
                        [num, den] = tfdata(K, 'v');
                        Kp = num(1);
                        Ki = num(2);
                        K = tf([Kp, Ki*0.01], [1, 0]);
                        
                    case 'PD'
                        % For PD, increase filtering and reduce high-frequency gain
                        [num, den] = tfdata(K, 'v');
                        K = tf(num*0.01, [epsilon*10, 1]);
                        
                    case 'PID'
                        % For PID, create a conservative controller using plant info
                        if plantInfo.isUnstable
                            p = plantInfo.poles;
                            unstable_poles = p(real(p) > 0);
                            
                            if ~isempty(unstable_poles)
                                [max_real, idx] = max(real(unstable_poles));
                                Kp = 2 * max_real;
                                Ki = 0.01 * Kp;
                                Kd = 5 * Kp / omega;
                                
                                K = tf([Kd, Kp, Ki], [epsilon*10*Kd, 1, 0]);
                            end
                        end
                end
                
                % Final stability check
                if isStateSpace
                    closed_loop_final = feedback(G_tf*K, 1);
                else
                    closed_loop_final = feedback(G*K, 1);
                end
                
                is_stable_final = all(real(pole(closed_loop_final)) < 0);
                
                if is_stable_final
                    details = [details, 'Successfully stabilized with structure-specific modifications.\n'];
                else
                    details = [details, 'WARNING: Could not stabilize system. Try a different design method.\n'];
                end
            end
        end
    catch ME
        details = [details, sprintf('\nError in closed-loop analysis: %s\n', ME.message)];
    end
    
    % Final controller properties
    try
        % Stability margins
        if isStateSpace
            [Gm, Pm, Wcg, Wcp] = margin(G_tf*K);
        else
            [Gm, Pm, Wcg, Wcp] = margin(G*K);
        end
        
        details = [details, sprintf('\nStability margins:\n')];
        details = [details, sprintf('  Gain margin: %.2f dB at %.4f rad/s\n', 20*log10(Gm), Wcg)];
        details = [details, sprintf('  Phase margin: %.2f deg at %.4f rad/s\n', Pm, Wcp)];
    catch
        details = [details, '\nCould not compute stability margins.\n'];
    end
    
    % Return final controller
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