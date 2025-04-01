function [K, details] = designPreStabilization(G, structure, options, plantInfo)
% DESIGNPRESTABILIZATION Controller design for highly unstable systems using pre-stabilization
%
% Implements a two-step controller design approach with enhanced state-space capabilities:
% 1. First stabilize the unstable plant with a simple controller
% 2. Then design a performance controller for the pre-stabilized system
%
% This method is particularly effective for highly unstable plants with
% multiple RHP poles, plants with both RHP poles and zeros, and plants with
% challenging dynamics.
%
% Inputs:
%   G        - Plant transfer function or state-space model
%   structure - Controller structure ('P', 'PI', 'PD', 'PID')
%   options   - Structure with design parameters:
%     .bandwidth  - Desired bandwidth in rad/s (default: 1)
%     .damping    - Desired damping ratio (default: 0.8)
%     .epsilon    - Filter parameter for D-term (default: 0.1)
%     .method     - Performance controller design method after pre-stabilization
%                  ('loop-shaping', 'pole-placement', 'h-infinity', default: 'loop-shaping')
%   plantInfo - Structure with plant analysis information
%
% Outputs:
%   K       - Designed controller as transfer function
%   details - Text description with design details

    % Start with detailed information about the method
    details = 'Pre-stabilization Controller Design Method\n';
    details = [details, '----------------------------------------\n'];
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
    
    if ~isfield(options, 'method')
        options.method = 'loop-shaping';
    end
    
    % Extract key parameters
    omega = options.bandwidth;
    zeta = options.damping;
    epsilon = options.epsilon;
    method = options.method;
    
    details = [details, sprintf('Desired bandwidth: %.4f rad/s\n', omega)];
    details = [details, sprintf('Desired damping ratio: %.4f\n', zeta)];
    details = [details, sprintf('Derivative filter coefficient: %.4f\n', epsilon)];
    details = [details, sprintf('Performance design method: %s\n', method)];
    
    % Detect if G is state-space model or transfer function
    isStateSpace = isa(G, 'ss');
    
    if isStateSpace
        details = [details, 'Using state-space model directly for controller design.\n'];
    end
    
    % Check if pre-stabilization is needed
    if ~plantInfo.isUnstable
        details = [details, 'Plant is already stable. Pre-stabilization not required.\n'];
        details = [details, 'Proceeding with direct controller design...\n'];
        
        % For stable plants, just use a standard design method
        switch lower(method)
            case 'loop-shaping'
                [K, method_details] = designLoopShaping(G, structure, options.phaseMargin, omega, epsilon, plantInfo);
            case 'pole-placement'
                [K, method_details] = designPolePlacement(G, structure, options, plantInfo);
            case 'h-infinity'
                robustness_opt = struct('robustness', 'Medium');
                [K, method_details] = designHInfinity(G, structure, robustness_opt, epsilon, plantInfo);
            otherwise
                [K, method_details] = designLoopShaping(G, structure, 45, omega, epsilon, plantInfo);
        end
        
        details = [details, method_details];
        return;
    end
    
    % Pre-stabilization phase
    details = [details, '\nSTEP 1: Pre-stabilization of unstable plant\n'];
    details = [details, '-------------------------------------\n'];
    
    % Get plant poles and zeros for analysis
    p = plantInfo.poles;
    try
        z = plantInfo.zeros;
    catch
        z = [];
    end
    
    % Identify unstable poles
    unstable_poles = p(real(p) > 0);
    details = [details, sprintf('Detected %d unstable poles:\n', length(unstable_poles))];
    
    for i = 1:length(unstable_poles)
        if imag(unstable_poles(i)) ~= 0
            details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(unstable_poles(i)), imag(unstable_poles(i)))];
        else
            details = [details, sprintf('  p%d = %.4f\n', i, real(unstable_poles(i)))];
        end
    end
    
    % Identify RHP zeros if any
    rhp_zeros = z(real(z) > 0);
    if ~isempty(rhp_zeros)
        details = [details, sprintf('Detected %d RHP zeros:\n', length(rhp_zeros))];
        
        for i = 1:length(rhp_zeros)
            if imag(rhp_zeros(i)) ~= 0
                details = [details, sprintf('  z%d = %.4f + %.4fi\n', i, real(rhp_zeros(i)), imag(rhp_zeros(i)))];
            else
                details = [details, sprintf('  z%d = %.4f\n', i, real(rhp_zeros(i)))];
            end
        end
    end
    
    % Initialize stabilizing controller
    K_stab = [];
    
    % Enhanced stabilization approach using state-space directly if available
    if isStateSpace
        details = [details, 'Using direct state-space approach for pre-stabilization.\n'];
        
        % Extract state-space matrices
        [A, B, C, D] = ssdata(G);
        n = size(A, 1);
        
        % Check controllability with robust numerical method
        Co = ctrb(A, B);
        sv = svd(Co);
        rank_tol = max(size(A)) * eps(norm(A));
        rank_Co = sum(sv > rank_tol);
        
        if rank_Co == n
            details = [details, 'System is fully controllable. Using state feedback for stabilization.\n'];
            
            % Design state feedback to place poles in LHP
            try
                % Create desired pole locations
                desired_poles = zeros(n, 1);
                p_idx = 1;
                
                % Handle unstable poles first
                for i = 1:length(unstable_poles)
                    pole_i = unstable_poles(i);
                    
                    if imag(pole_i) ~= 0
                        % Skip complex conjugate pairs (we'll handle them together)
                        if i < length(unstable_poles) && abs(pole_i - conj(unstable_poles(i+1))) < 1e-6
                            continue;
                        end
                        
                        if imag(pole_i) > 0 && p_idx <= n-1
                            % Complex pole - place a stable complex pair
                            magnitude = abs(pole_i);
                            desired_poles(p_idx) = -magnitude * 0.7 + 1j * imag(pole_i) * 0.7;
                            desired_poles(p_idx+1) = conj(desired_poles(p_idx));
                            p_idx = p_idx + 2;
                        end
                    else
                        % Real pole - reflect to LHP
                        if p_idx <= n
                            desired_poles(p_idx) = -abs(pole_i) * 1.5;
                            p_idx = p_idx + 1;
                        end
                    end
                end
                
                % Fill remaining poles
                for i = p_idx:n
                    desired_poles(i) = -omega * i;  % Additional poles based on bandwidth
                end
                
                details = [details, 'Desired pole locations for stabilizing feedback:\n'];
                for i = 1:length(desired_poles)
                    if imag(desired_poles(i)) ~= 0
                        details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(desired_poles(i)), imag(desired_poles(i)))];
                    else
                        details = [details, sprintf('  p%d = %.4f\n', i, real(desired_poles(i)))];
                    end
                end
                
                % Design state feedback gain
                try
                    F = place(A, B, desired_poles);
                    details = [details, 'Designed state feedback using place().\n'];
                catch
                    try
                        F = acker(A, B, desired_poles);
                        details = [details, 'Designed state feedback using acker().\n'];
                    catch
                        % Manual pole placement as fallback
                        F = -B' / (B*B') * (A - diag(desired_poles));
                        details = [details, 'Used manual pole placement as fallback.\n'];
                    end
                end
                
                % Design observer for state estimation
                obs_poles = desired_poles * 2.5;  % Faster observer poles
                
                try
                    L = place(A', C', obs_poles)';
                    details = [details, 'Designed observer using place().\n'];
                catch
                    try
                        L = acker(A', C', obs_poles)';
                        details = [details, 'Designed observer using acker().\n'];
                    catch
                        % Manual observer gain as fallback
                        L = -(A - diag(obs_poles))' / (C*C') * C';
                        details = [details, 'Used manual observer gain calculation as fallback.\n'];
                    end
                end
                
                % Create observer-based controller
                Ac = A - B*F - L*C;
                Bc = L;
                Cc = -F;
                Dc = 0;
                
                K_stab = ss(Ac, Bc, Cc, Dc);
                
                % Convert to transfer function if needed for later computation
                K_stab_tf = tf(K_stab);
                
                details = [details, 'Created observer-based stabilizing controller.\n'];
                
                % Verify stabilization
                try
                    G_stab = feedback(G, K_stab, 1, 1, -1);  % Negative feedback with state-space
                    stable_eigs = eig(G_stab.A);
                    is_stable = all(real(stable_eigs) < 0);
                    
                    if is_stable
                        details = [details, 'State-space design successfully stabilized the plant.\n'];
                    else
                        details = [details, 'WARNING: State-space design did not fully stabilize the plant.\n'];
                        
                        % Fallback to the transfer function approach
                        K_stab = [];
                    end
                catch ME
                    details = [details, sprintf('Error checking stability: %s\n', ME.message)];
                    
                    % Fallback to the transfer function approach
                    K_stab = [];
                end
            catch ME
                details = [details, sprintf('State feedback design failed: %s\n', ME.message)];
                
                % Fallback to transfer function approach
                K_stab = [];
            end
        else
            details = [details, sprintf('System not fully controllable (rank %d of %d).\n', rank_Co, n)];
            details = [details, 'Using transfer function approach for pre-stabilization.\n'];
        end
    end
    
    % If state-space approach failed or wasn't used, fall back to transfer function approach
    if isempty(K_stab)
        details = [details, 'Using transfer function approach for pre-stabilization.\n'];
        
        % Check for both RHP poles and zeros (most challenging case)
        if ~isempty(rhp_zeros)
            % For plants with both RHP poles and zeros, design is more challenging
            details = [details, 'Plant has both RHP poles and zeros - using advanced pre-stabilization.\n'];
            
            % Check for parity interlacing property (PIP) violations
            has_pip_violation = false;
            for i = 1:length(rhp_zeros)
                zero_i = rhp_zeros(i);
                poles_to_right = sum(real(unstable_poles) > real(zero_i));
                
                if poles_to_right > i
                    has_pip_violation = true;
                    details = [details, sprintf('PIP violation: %d poles to the right of zero at %.4f\n',poles_to_right, real(zero_i))];
                end
            end
            
            if has_pip_violation
                details = [details, 'Using specialized approach for PIP violations.\n'];
                
                % Create controller that addresses PIP violation
                % Approach: indirect stabilization with aggressive phase lead and filtering
                
                K_num = 1;
                K_den = 1;
                
                % First, add zeros to cancel unstable poles
                for i = 1:length(unstable_poles)
                    pole_i = unstable_poles(i);
                    
                    if imag(pole_i) ~= 0
                        % Skip complex conjugate pairs, we'll handle them together
                        if i < length(unstable_poles) && abs(pole_i - conj(unstable_poles(i+1))) < 1e-6
                            continue;
                        end
                        
                        % For complex poles, create zeros that cancel the pole pair
                        if imag(pole_i) > 0
                            real_part = real(pole_i);
                            imag_part = imag(pole_i);
                            
                            % Zero to cancel pole pair: s^2 - 2*re(p)*s + |p|^2
                            num_term = [1, -2*real_part, real_part^2 + imag_part^2];
                            
                            % Create stable poles in place of cancelled poles (further into LHP)
                            margin_factor = 3.0;  % Extra margin for PIP violation
                            den_term = [1, 2*margin_factor*abs(real_part), margin_factor^2*(real_part^2 + imag_part^2)];
                            
                            K_num = conv(K_num, num_term);
                            K_den = conv(K_den, den_term);
                        end
                    else
                        % For real unstable poles, add a zero to cancel and a stable pole
                        num_term = [1, -pole_i];
                        
                        margin_factor = 3.0;  % Extra margin for PIP violation
                        den_term = [1, margin_factor * abs(pole_i)];
                        
                        K_num = conv(K_num, num_term);
                        K_den = conv(K_den, den_term);
                    end
                end
                
                % Add aggressive filtering to handle RHP zeros
                for i = 1:length(rhp_zeros)
                    zero_i = rhp_zeros(i);
                    
                    if imag(zero_i) == 0
                        % For real RHP zeros, add low-pass filtering below the zero frequency
                        filter_freq = real(zero_i) / 3;  % Well below RHP zero
                        
                        filter_num = 1;
                        filter_den = [1, filter_freq];
                        
                        K_den = conv(K_den, filter_den);
                        
                        details = [details, sprintf('Added filtering at %.4f rad/s (below RHP zero at %.4f)\n',filter_freq, real(zero_i))];
                    end
                end
                
                % Calculate conservative gain
                max_real_part = max(real(unstable_poles));
                K_gain = max_real_part * 0.8;  % Conservative gain for PIP violations
                
                K_stab = tf(K_gain * K_num, K_den);
                details = [details, sprintf('Using gain factor of %.4f\n', K_gain)];
            else
                details = [details, 'No PIP violations detected. Using standard non-minimum phase stabilization.\n'];
                
                % For plants with RHP zeros but no PIP violations
                % We use phase lead networks carefully designed to respect RHP zero limitations
                
                K_num = 1;
                K_den = 1;
                
                % Create phase lead network
                min_rhp_zero = min(real(rhp_zeros));
                max_unstable_pole = max(real(unstable_poles));
                
                % Use phase lead centered between unstable pole and RHP zero
                if min_rhp_zero > max_unstable_pole
                    % If zero is to the right of all unstable poles, we have room to work
                    lead_zero = max_unstable_pole / 2;  % Zero halfway between origin and unstable pole
                    lead_pole = min_rhp_zero * 2;  % Pole well beyond RHP zero
                    
                    lead_num = [1, lead_zero];
                    lead_den = [1, lead_pole];
                    
                    K_num = conv(K_num, lead_num);
                    K_den = conv(K_den, lead_den);
                    
                    details = [details, sprintf('Added phase lead network with zero at %.4f and pole at %.4f\n',lead_zero, lead_pole)];
                else
                    % More challenging case - add multiple lead networks
                    % One for each unstable pole
                    
                    for i = 1:length(unstable_poles)
                        pole_i = unstable_poles(i);
                        
                        if imag(pole_i) == 0  % Handle real poles
                            % Phase lead network for each real unstable pole
                            lead_zero = pole_i * 0.8;  % Zero just to the left of unstable pole
                            lead_pole = -pole_i * 2;  % Stable pole at negative value
                            
                            lead_num = [1, -lead_zero];  % Note negative sign to place zero in RHP
                            lead_den = [1, -lead_pole];  % Negative sign to place pole in LHP
                            
                            K_num = conv(K_num, lead_num);
                            K_den = conv(K_den, lead_den);
                            
                            details = [details, sprintf('Added phase lead for pole at %.4f\n', pole_i)];
                        end
                    end
                    
                    % Add low-pass filtering to respect RHP zeros
                    filter_freq = min_rhp_zero / 3;
                    filter_den = [1, filter_freq];
                    K_den = conv(K_den, filter_den);
                    
                    details = [details, sprintf('Added filtering at %.4f rad/s (below RHP zero at %.4f)\n',filter_freq, min_rhp_zero)];
                end
                
                % Calculate appropriate gain
                max_real_part = max(real(unstable_poles));
                K_gain = max_real_part * 1.5;
                
                K_stab = tf(K_gain * K_num, K_den);
                details = [details, sprintf('Using gain factor of %.4f\n', K_gain)];
            end
        else
            % For plants with only unstable poles (no RHP zeros)
            details = [details, 'Plant has unstable poles only - using pole cancellation approach.\n'];
            
            % For multiple unstable poles, use pole-zero cancellation
            K_num = 1;
            K_den = 1;
            
            % Cancel each unstable pole with a zero
            for i = 1:length(unstable_poles)
                pole_i = unstable_poles(i);
                
                if imag(pole_i) ~= 0
                    % Skip complex conjugate pairs, we'll handle them together
                    if i < length(unstable_poles) && abs(pole_i - conj(unstable_poles(i+1))) < 1e-6
                        continue;
                    end
                    
                    % For complex poles, create zeros that cancel the pole pair
                    if imag(pole_i) > 0
                        real_part = real(pole_i);
                        imag_part = imag(pole_i);
                        
                        % Zero to cancel pole pair: s^2 - 2*re(p)*s + |p|^2
                        num_term = [1, -2*real_part, real_part^2 + imag_part^2];
                        
                        % Create stable poles in place of cancelled poles
                        margin_factor = 2.0;
                        den_term = [1, 2*margin_factor*abs(real_part), margin_factor^2*(real_part^2 + imag_part^2)];
                        
                        K_num = conv(K_num, num_term);
                        K_den = conv(K_den, den_term);
                    end
                else
                    % For real unstable poles, add a zero to cancel and a stable pole
                    num_term = [1, -pole_i];
                    
                    margin_factor = 2.0;
                    den_term = [1, margin_factor * abs(pole_i)];
                    
                    K_num = conv(K_num, num_term);
                    K_den = conv(K_den, den_term);
                end
            end
            
            % Add integral action for tracking if no integrator in plant
            if ~plantInfo.hasIntegrator && length(unstable_poles) <= 2  % Only for systems with few unstable poles
                % Add weak integral action
                int_freq = min(abs(real(unstable_poles))) / 10;
                int_den = [1, int_freq];
                K_den = conv(K_den, int_den);
                details = [details, sprintf('Added weak integral action with time constant %.4f\n', 1/int_freq)];
            end
            
            % Calculate appropriate gain
            max_real_part = max(real(unstable_poles));
            K_gain = max(1.0, max_real_part);
            
            K_stab = tf(K_gain * K_num, K_den);
            details = [details, sprintf('Using gain factor of %.4f\n', K_gain)];
        end
    end
    
    % Check if pre-stabilization is successful
    try
        % Handle different variable types
        if isStateSpace && isa(K_stab, 'ss')
            G_stab = feedback(G, K_stab, 1, 1, -1);  % Negative feedback with state-space
            stab_poles = eig(G_stab.A);
        else
            % If one or both are transfer functions, use standard feedback
            if isStateSpace
                G_tf = tf(G);  % Convert G to transfer function if needed
                G_stab = feedback(G_tf * K_stab, 1);
            else
                G_stab = feedback(G * K_stab, 1);
            end
            stab_poles = pole(G_stab);
        end
        
        is_prestab_successful = all(real(stab_poles) < 0);
        
        if is_prestab_successful
            details = [details, '\nPre-stabilization successful! Stabilized plant poles:\n'];
        else
            details = [details, '\nWARNING: Pre-stabilization not fully successful. Poles of inner loop:\n'];
            
            % Try adjusting the controller gain
            if ~isStateSpace || ~isa(K_stab, 'ss')
                % Get controller data
                [num, den] = tfdata(K_stab, 'v');
                
                % Try different gain scaling factors
                scale_factors = [0.5, 0.2, 0.1, 0.05, 0.01];
                
                for factor = scale_factors
                    K_test = tf(num * factor, den);
                    
                    if isStateSpace
                        G_tf = tf(G);
                        G_test = feedback(G_tf * K_test, 1);
                    else
                        G_test = feedback(G * K_test, 1);
                    end
                    
                    if all(real(pole(G_test)) < 0)
                        K_stab = K_test;
                        details = [details, sprintf('System stabilized by scaling controller gain by %.2f\n', factor)];
                        
                        if isStateSpace
                            G_stab = feedback(G_tf * K_stab, 1);
                        else
                            G_stab = feedback(G * K_stab, 1);
                        end
                        
                        stab_poles = pole(G_stab);
                        is_prestab_successful = true;
                        break;
                    end
                end
            end
        end
        
        for i = 1:length(stab_poles)
            if imag(stab_poles(i)) ~= 0
                details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(stab_poles(i)), imag(stab_poles(i)))];
            else
                details = [details, sprintf('  p%d = %.4f\n', i, real(stab_poles(i)))];
            end
        end
        
        % If pre-stabilization failed, try alternative approach
        if ~is_prestab_successful
            details = [details, 'Initial pre-stabilization unsuccessful. Trying alternative approach...\n'];
            
            % Try more aggressive stabilization techniques
            try
                % Create a very simple, extremely conservative controller
                if ~isempty(rhp_zeros)
                    % For systems with RHP zeros, build a controller respecting these zeros
                    min_rhp_zero = min(real(rhp_zeros));
                    rolloff_freq = min_rhp_zero / 5;  % Roll off well below RHP zero
                    
                    % Simple lead controller with heavy filtering
                    K_stab = tf([1, 0.1], [1/rolloff_freq, 1]);
                    details = [details, sprintf('Created conservative lead controller with rolloff at %.4f rad/s\n', rolloff_freq)];
                else
                    % For systems without RHP zeros, use a simple lead controller
                    K_stab = tf([1, 0.1], [0.01, 1]);
                    details = [details, 'Created conservative lead controller with high gain at high frequencies\n'];
                end
                
                % Try different gain levels for stabilization
                [num, den] = tfdata(K_stab, 'v');
                is_prestab_successful = false;
                
                for gain = [5, 2, 1, 0.5, 0.1, 0.05, 0.01, 0.005]
                    K_test = tf(num * gain, den);
                    
                    if isStateSpace
                        G_tf = tf(G);
                        G_test = feedback(G_tf * K_test, 1);
                    else
                        G_test = feedback(G * K_test, 1);
                    end
                    
                    if all(real(pole(G_test)) < 0)
                        K_stab = K_test;
                        details = [details, sprintf('Found stabilizing controller with gain = %.4f\n', gain)];
                        
                        % Update G_stab with the successful controller
                        if isStateSpace
                            G_stab = feedback(G_tf * K_stab, 1);
                        else
                            G_stab = feedback(G * K_stab, 1);
                        end
                        
                        stab_poles = pole(G_stab);
                        is_prestab_successful = true;
                        break;
                    end
                end
            catch ME
                details = [details, sprintf('Alternative approach failed: %s\n', ME.message)];
                details = [details, 'Proceeding with best effort stabilizing controller.\n'];
                
                % Use the last K_stab we tried, even if not fully successful
            end
        end
        
        % Continue with the second step only if stabilization was successful
        if is_prestab_successful
            details = [details, '\nSTEP 2: Design performance controller for pre-stabilized plant\n'];
            details = [details, '--------------------------------------------------------\n'];
            
            % Now design a performance controller for the pre-stabilized plant
            % Create a new plantInfo structure for the stabilized plant
            plantInfo_stab = analyzePlant(G_stab);
            
            % Adjust bandwidth if needed based on pre-stabilized dynamics
            omega_adj = min(omega, getBandwidth(G_stab, plantInfo_stab) * 0.5);
            if omega_adj < omega
                details = [details, sprintf('Adjusted bandwidth to %.4f rad/s based on pre-stabilized plant dynamics.\n', omega_adj)];
            end
            
            % Design performance controller using selected method
            switch lower(method)
                case 'loop-shaping'
                    [K_perf, loop_details] = designLoopShaping(G_stab, structure, 45, omega_adj, epsilon, plantInfo_stab);
                    details = [details, loop_details];
                    
                case 'pole-placement'
                    options_adj = options;
                    options_adj.bandwidth = omega_adj;
                    [K_perf, pp_details] = designPolePlacement(G_stab, structure, options_adj, plantInfo_stab);
                    details = [details, pp_details];
                    
                case 'h-infinity'
                    robustness_opt = struct('robustness', 'Medium');
                    [K_perf, hinf_details] = designHInfinity(G_stab, structure, robustness_opt, epsilon, plantInfo_stab);
                    details = [details, hinf_details];
                    
                otherwise
                    [K_perf, default_details] = designLoopShaping(G_stab, structure, 45, omega_adj, epsilon, plantInfo_stab);
                    details = [details, default_details];
            end
            
            % Combine pre-stabilizing and performance controllers
            if isa(K_stab, 'ss') && ~isa(K_perf, 'ss')
                % Convert K_stab to transfer function for easier combination
                K_stab_tf = tf(K_stab);
                K = series(K_perf, K_stab_tf);
            elseif ~isa(K_stab, 'ss') && isa(K_perf, 'ss')
                % Convert K_perf to transfer function
                K_perf_tf = tf(K_perf);
                K = series(K_perf_tf, K_stab);
            else
                % Both same type - either both ss or both tf
                K = series(K_perf, K_stab);
            end
            
            % Simplify the combined controller if possible
            try
                K_simple = minreal(K, 0.01);
                details = [details, sprintf('\nReduced controller order from %d to %d through minimization.\n',order(K), order(K_simple))];
                K = K_simple;
            catch
                details = [details, '\nCould not simplify combined controller.\n'];
            end
            
            % Analyze final closed-loop system
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
                    details = [details, '\nFinal closed-loop system is stable! Actual poles:\n'];
                else
                    details = [details, '\nWARNING: Final closed-loop system is UNSTABLE! Actual poles:\n'];
                end
                
                for i = 1:length(cl_poles)
                    if imag(cl_poles(i)) ~= 0
                        details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(cl_poles(i)), imag(cl_poles(i)))];
                    else
                        details = [details, sprintf('  p%d = %.4f\n', i, real(cl_poles(i)))];
                    end
                end
                
                % If final system is unstable, make adjustments
                if ~is_stable
                    details = [details, '\nAttempting to stabilize final controller by gain adjustment...\n'];
                    
                    % Try reducing the gain
                    scales = [0.5, 0.2, 0.1, 0.05, 0.01];
                    for scale = scales
                        K_scaled = K * scale;
                        
                        if isStateSpace
                            closed_loop_scaled = feedback(G_tf * K_scaled, 1);
                        else
                            closed_loop_scaled = feedback(G * K_scaled, 1);
                        end
                        
                        if all(real(pole(closed_loop_scaled)) < 0)
                            K = K_scaled;
                            details = [details, sprintf('Final system stabilized by scaling controller gain by %.2f\n', scale)];
                            break;
                        end
                    end
                    
                    % Final check
                    if isStateSpace
                        final_cl = feedback(G_tf*K, 1);
                    else
                        final_cl = feedback(G*K, 1);
                    end
                    
                    final_stable = all(real(pole(final_cl)) < 0);
                    
                    if ~final_stable
                        details = [details, 'WARNING: Could not stabilize final system. Using pre-stabilizing controller only.\n'];
                        K = K_stab;
                    end
                end
                
                % Calculate stability margins for final controller
                try
                    if isStateSpace
                        [Gm, Pm, Wcg, Wcp] = margin(G_tf*K);
                    else
                        [Gm, Pm, Wcg, Wcp] = margin(G*K);
                    end
                    
                    details = [details, sprintf('\nStability margins of final controller:\n')];
                    details = [details, sprintf('  Gain margin: %.2f dB at %.4f rad/s\n', 20*log10(Gm), Wcg)];
                    details = [details, sprintf('  Phase margin: %.2f deg at %.4f rad/s\n', Pm, Wcp)];
                catch
                    details = [details, '\nCould not compute stability margins for final controller.\n'];
                end
            catch ME
                details = [details, sprintf('\nError in final closed-loop analysis: %s\n', ME.message)];
            end
        else
            % If all pre-stabilization attempts fail, return a simplified controller
            details = [details, 'WARNING: All pre-stabilization attempts failed.\n'];
            details = [details, 'Providing simple structure-specific controller. Manual tuning strongly recommended.\n'];
            
            % Create a basic controller based on structure
            switch structure
                case 'P'
                    K = tf(0.5, 1);
                case 'PI'
                    K = tf([1, 0.1], [1, 0]);
                case 'PD'
                    K = tf([0.5, 0.1], [0.01, 1]);
                case 'PID'
                    K = tf([0.5, 1, 0.1], [0.01, 1, 0]);
                otherwise
                    K = tf(0.1, 1);
            end
        end
    catch ME
        details = [details, sprintf('Error during pre-stabilization: %s\n', ME.message)];
        
        % Emergency fallback controller based on structure
        switch structure
            case 'P'
                K = tf(0.1, 1);
            case 'PI'
                K = tf([0.1, 0.01], [1, 0]);
            case 'PD'
                K = tf([0.1, 0.01], [0.01, 1]);
            case 'PID'
                K = tf([0.1, 0.1, 0.01], [0.01, 1, 0]);
            otherwise
                K = tf(0.1, 1);
        end
        
        details = [details, 'Emergency fallback controller created. Manual tuning required.\n'];
    end
    
    return;
end

% Function to get bandwidth estimate
function w_bandwidth = getBandwidth(G, plantInfo)
    % Try to estimate system bandwidth
    try
        w_bandwidth = bandwidth(G);
        if ~isnan(w_bandwidth)
            return;
        end
    catch
        % Continue to alternative methods
    end
    
    % Alternative method: Find -3dB point from Bode plot
    try
        w = logspace(-3, 3, 200);
        [mag, ~] = bode(G, w);
        mag = squeeze(mag);
        
        % Find DC gain or low-frequency gain
        try
            dc_gain = dcgain(G);
        catch
            dc_gain = mag(1);
        end
        
        % Find -3dB frequency
        threshold = 0.7071 * abs(dc_gain);  % -3dB = 0.7071x
        idx = find(mag < threshold, 1);
        
        if ~isempty(idx) && idx > 1
            w_bandwidth = w(idx);
            return;
        end
    catch
        % Continue to next method
    end
    
    % Use dominant poles
    p = plantInfo.poles;
    stable_poles = p(real(p) < 0);
    
    if ~isempty(stable_poles)
        % Find the slowest stable pole
        [~, idx] = min(abs(real(stable_poles)));
        w_bandwidth = abs(real(stable_poles(idx))) * 2;  % Rule of thumb
        return;
    end
    
    % Default fallback
    w_bandwidth = 1.0;
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