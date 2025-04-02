function [K, details] = designEnhancedStateFeedback(G, structure, options, plantInfo)
% DESIGNENHANCEDSTATEFEEDBACK Advanced controller design using state feedback with observer
%
% Implements robust state feedback with observer specifically designed for
% unstable and challenging systems. The method combines pole placement with
% state estimation and includes special handling for high-order unstable systems.
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
        % Convert to state space if not already
        try
            % Try the direct conversion first
            G_ss = ss(G);
            isStateSpace = true;
            details = [details, 'Successfully converted to state-space representation.\n'];
        catch ME
            % If direct conversion fails, try a different approach for problematic cases
            details = [details, sprintf('Direct state-space conversion failed: %s\n', ME.message)];
            details = [details, 'Attempting alternative conversion approach...\n'];
            
            try
                % If there's a pole at the origin (integrator), handle carefully
                if any(abs(pole(G)) < 1e-6)
                    details = [details, 'Plant contains integrator(s). Using modified approach.\n'];
                    
                    % For transfer functions with integrators, first cancel common factors
                    G_clean = minreal(G);
                    
                    % Then, manually check for s in numerator and denominator
                    [num, den] = tfdata(G_clean, 'v');
                    
                    % If denominator has a zero at the end (s factor), remove it
                    if abs(den(end)) < 1e-10
                        den = den(1:end-1);
                    end
                    
                    % If both numerator and denominator have s factors, remove one from each
                    if abs(num(end)) < 1e-10 && abs(den(end)) < 1e-10
                        num = num(1:end-1);
                        den = den(1:end-1);
                    end
                    
                    % Create a new transfer function without the common s factor
                    G_reduced = tf(num, den);
                    details = [details, 'Manually removed common s factor.\n'];
                    
                    % Try converting the reduced system
                    G_ss = ss(G_reduced);
                else
                    % Try to balance the system to improve numerical conditioning
                    [num, den] = tfdata(G, 'v');
                    sys_tf = tf(num, den);
                    G_ss = balreal(ss(sys_tf));
                end
                
                isStateSpace = true;
            catch ME2
                details = [details, sprintf('Alternative conversion failed: %s\n', ME2.message)];
                
                % Last resort: try to model the system numerically
                details = [details, 'Attempting numerical approximation...\n'];
                
                try
                    % Get zeros and poles
                    z = zero(G);
                    p = pole(G);
                    k = 1; % Use unity gain if dcgain fails
                    
                    try
                        k = dcgain(G);
                        if isinf(k) || isnan(k)
                            k = 1;
                        end
                    catch
                        % Keep k = 1
                    end
                    
                    % Create a balanced representation by hand
                    A = diag(p);
                    B = ones(length(p), 1);
                    C = ones(1, length(p));
                    D = 0;
                    
                    % Try to improve the model if we have zeros
                    if ~isempty(z)
                        % Try to incorporate zeros but keep it simple
                        C = [ones(1, length(p)-length(z)), zeros(1, length(z))];
                    end
                    
                    % Use balanced state-space model
                    G_ss = ss(A, B, C, D);
                    isStateSpace = true;
                    
                    details = [details, 'Using numerically balanced model.\n'];
                catch ME3
                    details = [details, sprintf('All conversion attempts failed: %s\n', ME3.message)];
                    details = [details, 'Using alternative controller design method.\n'];
                    
                    % Fallback to using a completely different method
                    [K, fallback_details] = designPreStabilization(G, structure, options, plantInfo);
                    details = [details, 'CONVERSION FAILED: Using Pre-stabilization method instead.\n\n'];
                    details = [details, fallback_details];
                    return;
                end
            end
        end
    end
    
    % Extract state-space matrices
    A = G_ss.A;
    B = G_ss.B;
    C = G_ss.C;
    D = G_ss.D;
    
    % Get system dimensions
    n = size(A, 1); % Number of states
    
    % Check controllability and observability using SVD for numerical robustness
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
        A_reg = A + 1e-6 * eye(n);
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
        
        % Scale bandwidth for highly unstable systems
        if ~isempty(unstable_idx)
            max_real_part = max(real(p(unstable_idx)));
            
            if max_real_part > 5 || num_unstable > 1
                % For highly unstable systems, use very conservative bandwidth
                conservative_omega = max(omega, max_real_part/2);
                details = [details, sprintf('System is highly unstable. Using conservative bandwidth %.4f rad/s.\n', conservative_omega)];
            else
                conservative_omega = max(omega, 2 * max_real_part);
                details = [details, sprintf('Adjusted bandwidth to %.4f rad/s based on unstable poles.\n', conservative_omega)];
            end
        else
            conservative_omega = omega;
        end
        
        % Place poles in the LHP with appropriate damping
        desired_poles = zeros(n, 1);
        
        % Place unstable poles in LHP with sufficient damping
        for i = 1:n
            if i <= length(unstable_idx)
                idx = unstable_idx(i);
                
                if imag(p(idx)) ~= 0
                    % For complex unstable poles, maintain frequency but add damping
                    freq = abs(p(idx));
                    desired_poles(i) = -zeta * freq + 1j * freq * sqrt(1 - zeta^2);
                    
                    % Make sure we handle the conjugate pair
                    if i < n && imag(p(idx)) > 0
                        desired_poles(i+1) = conj(desired_poles(i));
                    end
                else
                    % For real unstable poles, move to LHP with appropriate spacing
                    desired_poles(i) = -conservative_omega * (1 + 0.2 * (i-1));
                end
            else
                % For remaining states, place poles at multiples of bandwidth
                desired_poles(i) = -conservative_omega * (1.5 + 0.3 * (i-length(unstable_idx)));
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
            sf_poles = dominant_poles(1:n);
        end
        
        desired_poles = sf_poles';
    end
    
    details = [details, 'Desired state feedback poles:\n'];
    for i = 1:length(desired_poles)
        if imag(desired_poles(i)) ~= 0
            details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(desired_poles(i)), imag(desired_poles(i)))];
        else
            details = [details, sprintf('  p%d = %.4f\n', i, real(desired_poles(i)))];
        end
    end
    
    % Calculate observer poles - typically faster than closed-loop poles
    % Apply robustness setting to observer poles
    switch robustness
        case 'Low'
            obs_factor = 2.0;  % Faster observer, less filtering
        case 'High'
            obs_factor = 4.0;  % Slower observer, more filtering
        otherwise
            obs_factor = 2.5;  % Medium setting
    end
    
    % For highly unstable systems, use even faster observer
    if plantInfo.isUnstable
        p = eig(A);
        unstable_idx = find(real(p) > 0);
        
        if ~isempty(unstable_idx)
            max_real_part = max(real(p(unstable_idx)));
            
            if max_real_part > 5 || length(unstable_idx) > 1
                obs_factor = obs_factor * 1.5;
                details = [details, 'System is highly unstable. Using faster observer dynamics.\n'];
            end
        end
    end
    
    obs_poles = desired_poles * obs_factor;
    
    details = [details, sprintf('\nObserver poles (placed at a factor of %.1f times control poles):\n', obs_factor)];
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
        K_sf = place(A, B, desired_poles);
        details = [details, 'Successfully computed state feedback gain K using place().\n'];
    catch ME
        details = [details, sprintf('Standard place() failed: %s\n', ME.message)];
        
        try
            % Try acker() as an alternative
            K_sf = acker(A, B, desired_poles);
            details = [details, 'Computed state feedback gain K using acker().\n'];
        catch ME2
            details = [details, sprintf('Acker method failed: %s\n', ME2.message)];
            
            % Use a more robust manual pole placement approach
            details = [details, 'Using robust manual pole placement.\n'];
            
            % Manual pole placement based on SVD for better numerical properties
            K_sf = zeros(1, n);
            
            % For unstable systems, focus on moving unstable poles first
            if plantInfo.isUnstable
                p = eig(A);
                unstable_idx = find(real(p) > 0);
                
                for i = 1:length(unstable_idx)
                    idx = unstable_idx(i);
                    
                    % Create a shift vector focused on this eigenvector
                    [V, ~] = eig(A);
                    shift_vector = V(:, idx);
                    
                    % Calculate how much to shift (to desired pole location)
                    new_pole = -abs(real(p(idx))) * 2;
                    shift_amount = new_pole - real(p(idx));
                    
                    % Update the gain
                    K_add = shift_amount * shift_vector' * pinv(B);
                    K_sf = K_sf + K_add;
                end
                
                details = [details, 'Created stabilizing gain focusing on unstable poles.\n'];
            else
                % For stable systems, use a simpler approach
                for i = 1:n
                    K_sf = K_sf + desired_poles(i) * pinv(B) * pinv(A - desired_poles(i)*eye(n));
                end
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
            
            % Use a more robust manual observer gain calculation
            details = [details, 'Using robust manual observer gain calculation.\n'];
            
            % Calculate observer gain using SVD for numerical robustness
            L = zeros(n, 1);
            
            if plantInfo.isUnstable
                % For unstable plants, ensure observer is fast enough
                p = eig(A);
                unstable_idx = find(real(p) > 0);
                
                for i = 1:length(unstable_idx)
                    idx = unstable_idx(i);
                    
                    % Create a shift vector focused on this eigenvector
                    [V, ~] = eig(A');
                    shift_vector = V(:, idx);
                    
                    % Calculate how much to shift (to desired pole location)
                    new_pole = -abs(real(p(idx))) * 4;
                    shift_amount = new_pole - real(p(idx));
                    
                    % Update the gain
                    L_add = shift_amount * pinv(C') * shift_vector;
                    L = L + L_add;
                end
            else
                % For stable systems, use a simpler approach
                for i = 1:n
                    L = L + obs_poles(i) * pinv(C') * pinv(A' - obs_poles(i)*eye(n));
                end
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
        details = [details, sprintf('Reduced controller order from %d to %d through minimization.\n', order(controller_tf), order(simplified_tf))];
        controller_tf = simplified_tf;
    catch
        details = [details, 'Could not simplify controller transfer function.\n'];
    end
    
    % Create the simplified controller in the requested structure
    try
        % Use frequency response analysis for robust parameter extraction
        w = logspace(-3, log10(omega*10), 200);
        [mag, phase] = bode(controller_tf, w);
        mag = squeeze(mag);
        phase = squeeze(phase);
        
        switch structure
            case 'P'
                % Extract proportional term at crossover frequency
                idx_mid = ceil(length(w)/2);
                Kp = mag(idx_mid);
                
                % Limit to reasonable values
                Kp = min(max(abs(Kp), 0.1), 100);
                
                K = tf(Kp, 1);
                details = [details, sprintf('Created P controller: Kp = %.4f\n', Kp)];
                
            case 'PI'
                % Check for integral action in the frequency response
                has_integral = (phase(1) < -45);
                
                if has_integral
                    details = [details, 'Detected integral action in observer-based controller.\n'];
                    
                    % Find frequency where phase crosses -45 degrees
                    idx_45 = find(phase > -45, 1);
                    if isempty(idx_45)
                        idx_45 = 2;
                    end
                    w_i = w(idx_45);
                    
                    % Get gain at mid-frequency for proportional term
                    idx_mid = ceil(length(w)/2);
                    Kp = mag(idx_mid);
                    
                    % Calculate integral gain based on phase crossover
                    Ki = Kp * w_i / 5;
                else
                    details = [details, 'No clear integral action detected. Adding appropriate integral term.\n'];
                    
                    % Extract proportional gain at mid-frequency
                    idx_mid = ceil(length(w)/2);
                    Kp = mag(idx_mid);
                    
                    % Add integral action based on bandwidth
                    Ki = Kp * omega / 10;
                    
                    if plantInfo.hasIntegrator
                        Ki = Ki * 0.5;  % Reduce integral action for plants with integrator
                    end
                end
                
                % For unstable systems, reduce integral action
                if plantInfo.isUnstable
                    p = eig(A);
                    unstable_idx = find(real(p) > 0);
                    
                    if ~isempty(unstable_idx)
                        max_real_part = max(real(p(unstable_idx)));
                        
                        if max_real_part > 5 || length(unstable_idx) > 1
                            Ki = Ki * 0.1;  % Greatly reduce integral action for highly unstable plants
                            details = [details, 'Greatly reduced integral action for highly unstable plant.\n'];
                        else
                            Ki = Ki * 0.5;  % Moderately reduce integral action for unstable plants
                            details = [details, 'Reduced integral action for unstable plant.\n'];
                        end
                    end
                end
                
                % Limit to reasonable values
                Kp = min(max(abs(Kp), 0.1), 100);
                Ki = min(max(abs(Ki), 0.01), 50);
                
                K = tf([Kp, Ki], [1, 0]);
                details = [details, sprintf('Created PI controller: Kp = %.4f, Ki = %.4f\n', Kp, Ki)];
                
            case 'PD'
                % Check for derivative action in the frequency response
                has_derivative = any(phase > 10);
                
                if has_derivative
                    details = [details, 'Detected derivative action in observer-based controller.\n'];
                    
                    % Find frequency of maximum phase (derivative action)
                    [~, idx_max] = max(phase);
                    w_d = w(idx_max);
                    
                    % Get proportional gain at mid-frequency
                    idx_mid = ceil(length(w)/2);
                    Kp = mag(idx_mid);
                    
                    % Calculate derivative gain
                    Kd = Kp / w_d;
                else
                    details = [details, 'No clear derivative action detected. Adding appropriate derivative term.\n'];
                    
                    % Extract proportional gain at mid-frequency
                    idx_mid = ceil(length(w)/2);
                    Kp = mag(idx_mid);
                    
                    % Add derivative action based on bandwidth
                    Kd = Kp / omega;
                    
                    if plantInfo.hasRHPZeros
                        Kd = Kd * 0.5;  % Reduce derivative action for non-minimum phase plants
                    end
                end
                
                % Limit to reasonable values
                Kp = min(max(abs(Kp), 0.1), 100);
                Kd = min(max(abs(Kd), 0.01), 50);
                
                % Add filtering for derivative term
                Td = Kd / Kp;
                K = tf([Kd, Kp], [epsilon*Td, 1]);
                details = [details, sprintf('Created PD controller: Kp = %.4f, Kd = %.4f\n', Kp, Kd)];
                
            case 'PID'
                % Check for both integral and derivative actions
                has_integral = (phase(1) < -45);
                has_derivative = any(phase > 10);
                
                if has_integral && has_derivative
                    details = [details, 'Detected both integral and derivative actions in controller.\n'];
                    
                    % Find frequency where phase crosses -45 degrees (integral action)
                    idx_45 = find(phase > -45, 1);
                    if isempty(idx_45)
                        idx_45 = 2;
                    end
                    w_i = w(idx_45);
                    
                    % Find frequency of maximum phase (derivative action)
                    [~, idx_max] = max(phase);
                    w_d = w(idx_max);
                    
                    % Get proportional gain at mid-frequency between integral and derivative
                    idx_mid = ceil((idx_45 + idx_max) / 2);
                    Kp = mag(idx_mid);
                    
                    % Calculate integral and derivative gains
                    Ki = Kp * w_i / 5;
                    Kd = Kp / w_d;
                else
                    details = [details, 'Incomplete PID behavior detected. Creating balanced PID controller.\n'];
                    
                    % Extract proportional gain at mid-frequency
                    idx_mid = ceil(length(w)/2);
                    Kp = mag(idx_mid);
                    
                    % Add integral and derivative actions based on bandwidth
                    Ki = Kp * omega / 10;
                    Kd = Kp / omega;
                    
                    if plantInfo.hasIntegrator
                        Ki = Ki * 0.5;  % Reduce integral action for plants with integrator
                    end
                    
                    if plantInfo.hasRHPZeros
                        Kd = Kd * 0.5;  % Reduce derivative action for non-minimum phase plants
                    end
                end
                
                % For unstable systems, adjust gains further
                if plantInfo.isUnstable
                    p = eig(A);
                    unstable_idx = find(real(p) > 0);
                    
                    if ~isempty(unstable_idx)
                        max_real_part = max(real(p(unstable_idx)));
                        
                        if max_real_part > 5 || length(unstable_idx) > 1
                            Ki = Ki * 0.1;  % Greatly reduce integral action for highly unstable plants
                            details = [details, 'Greatly reduced integral action for highly unstable plant.\n'];
                        else
                            Ki = Ki * 0.5;  % Moderately reduce integral action for unstable plants
                            details = [details, 'Reduced integral action for unstable plant.\n'];
                        end
                    end
                end
                
                % Limit to reasonable values
                Kp = min(max(abs(Kp), 0.1), 100);
                Ki = min(max(abs(Ki), 0.01), 50);
                Kd = min(max(abs(Kd), 0.01), 50);
                
                K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                details = [details, sprintf('Created PID controller: Kp = %.4f, Ki = %.4f, Kd = %.4f\n', Kp, Ki, Kd)];
                
            otherwise
                error('Unsupported controller structure');
        end
    catch ME
        details = [details, sprintf('Error extracting controller parameters: %s\n', ME.message)];
        details = [details, 'Using fallback approach for parameter extraction.\n'];
        
        % Fallback to using DC gain and order for extraction
        try
            Kp = abs(dcgain(controller_tf));
            if isnan(Kp) || isinf(Kp)
                Kp = 1.0;  % Default value
            end
            
            % Limit gain to reasonable values
            Kp = min(max(Kp, 0.1), 100);
            
            % Create controller based on structure
            switch structure
                case 'P'
                    K = tf(Kp, 1);
                    details = [details, sprintf('Created basic P controller: Kp = %.4f\n', Kp)];
                    
                case 'PI'
                    Ki = Kp * omega / 10;
                    
                    if plantInfo.hasIntegrator
                        Ki = Ki * 0.5;
                    end
                    
                    if plantInfo.isUnstable
                        Ki = Ki * 0.2;
                    end
                    
                    K = tf([Kp, Ki], [1, 0]);
                    details = [details, sprintf('Created basic PI controller: Kp = %.4f, Ki = %.4f\n', Kp, Ki)];
                    
                case 'PD'
                    Kd = Kp / omega;
                    
                    if plantInfo.hasRHPZeros
                        Kd = Kd * 0.5;
                    end
                    
                    Td = Kd / Kp;
                    K = tf([Kd, Kp], [epsilon*Td, 1]);
                    details = [details, sprintf('Created basic PD controller: Kp = %.4f, Kd = %.4f\n', Kp, Kd)];
                    
                case 'PID'
                    Ki = Kp * omega / 10;
                    Kd = Kp / omega;
                    
                    if plantInfo.hasIntegrator
                        Ki = Ki * 0.5;
                    end
                    
                    if plantInfo.hasRHPZeros
                        Kd = Kd * 0.5;
                    end
                    
                    if plantInfo.isUnstable
                        Ki = Ki * 0.2;
                    end
                    
                    K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                    details = [details, sprintf('Created basic PID controller: Kp = %.4f, Ki = %.4f, Kd = %.4f\n', Kp, Ki, Kd)];
                    
                otherwise
                    K = tf(Kp, 1);
                    details = [details, 'Using fallback P controller.\n'];
            end
        catch ME2
            details = [details, sprintf('Fallback extraction also failed: %s\n', ME2.message)];
            
            % Create a default controller if all else fails
            switch structure
                case 'P'
                    K = tf(1.0, 1);
                case 'PI'
                    K = tf([1.0, 0.1], [1, 0]);
                case 'PD'
                    K = tf([0.1, 1.0], [epsilon*0.1, 1]);
                case 'PID'
                    K = tf([0.1, 1.0, 0.1], [epsilon*0.1, 1, 0]);
                otherwise
                    K = tf(1.0, 1);
            end
            
            details = [details, 'Using default controller parameters.\n'];
        end
    end
    
    % Verify closed-loop stability with the extracted controller
    try
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
        
        % If unstable, try gain reduction
        if ~is_stable
            details = [details, '\nAttempting to stabilize by gain reduction...\n'];
            
            [num, den] = tfdata(K, 'v');
            stabilized = false;
            
            % Try various scaling factors
            for scale = [0.5, 0.2, 0.1, 0.05, 0.01]
                K_scaled = tf(num * scale, den);
                closed_loop_scaled = feedback(G * K_scaled, 1);
                
                if all(real(pole(closed_loop_scaled)) < 0)
                    K = K_scaled;
                    stabilized = true;
                    details = [details, sprintf('System stabilized with gain scaling factor: %.3f\n', scale)];
                    
                    % Show stabilized poles
                    cl_poles = pole(closed_loop_scaled);
                    details = [details, 'Stabilized closed-loop poles:\n'];
                    for i = 1:length(cl_poles)
                        if imag(cl_poles(i)) ~= 0
                            details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(cl_poles(i)), imag(cl_poles(i)))];
                        else
                            details = [details, sprintf('  p%d = %.4f\n', i, real(cl_poles(i)))];
                        end
                    end
                    
                    break;
                end
            end
            
            if ~stabilized
                details = [details, 'WARNING: Could not stabilize by gain scaling.\n'];
                
                % Try structure-specific modifications
                details = [details, 'Attempting structure-specific modifications...\n'];
                
                switch structure
                    case 'P'
                        % For P, just try an extremely low gain
                        K = tf(0.01, 1);
                        
                    case 'PI'
                        % For PI, drastically reduce integral action
                        [num, den] = tfdata(K, 'v');
                        Kp = num(1) * 0.1;
                        Ki = num(2) * 0.01;
                        K = tf([Kp, Ki], [1, 0]);
                        
                    case 'PD'
                        % For PD, increase filtering and reduce derivative action
                        [num, den] = tfdata(K, 'v');
                        Kp = num(2) * 0.1;
                        Kd = num(1) * 0.05;
                        Td = Kd / Kp;
                        K = tf([Kd, Kp], [epsilon*10*Td, 1]);
                        
                    case 'PID'
                        % For PID, create a conservative controller
                        if plantInfo.isUnstable
                            % For unstable plants, focus on stabilization
                            Kp = 0.1;
                            Ki = 0.001;
                            Kd = 0.5;
                            K = tf([Kd, Kp, Ki], [epsilon*10*Kd, 1, 0]);
                        else
                            % For stable plants, use a balanced approach
                            Kp = 0.1;
                            Ki = 0.01;
                            Kd = 0.1;
                            K = tf([Kd, Kp, Ki], [epsilon*5*Kd, 1, 0]);
                        end
                end
                
                % Check if modifications helped
                closed_loop_final = feedback(G*K, 1);
                is_stable_final = all(real(pole(closed_loop_final)) < 0);
                
                if is_stable_final
                    details = [details, 'Successfully stabilized with structure-specific modifications.\n'];
                else
                    details = [details, 'WARNING: All stabilization attempts failed. This system is extremely challenging.\n'];
                    
                    % Try pre-stabilization as a last resort
                    details = [details, 'Attempting pre-stabilization as last resort...\n'];
                    
                    try
                        [K_prestab, ~] = preStabilize(G, plantInfo);
                        
                        closed_loop_prestab = feedback(G*K_prestab, 1);
                        if all(real(pole(closed_loop_prestab)) < 0)
                            K = K_prestab;
                            details = [details, 'Using pre-stabilization controller as final fallback.\n'];
                        else
                            details = [details, 'Pre-stabilization also failed. Manual tuning required.\n'];
                        end
                    catch
                        details = [details, 'Pre-stabilization fallback failed. Manual tuning required.\n'];
                    end
                end
            end
        end
    catch ME
        details = [details, sprintf('\nError in closed-loop analysis: %s\n', ME.message)];
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