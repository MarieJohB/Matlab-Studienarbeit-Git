function [K, details] = designEmergencyController(G, structure, options, plantInfo)
    % DESIGNEMERGENCYCONTROLLER Design a robust emergency controller for difficult plants
    % This specialized controller is designed to stabilize extremely challenging
    % control plants when all other methods fail
    %
    % Inputs:
    %   G         - Plant transfer function
    %   structure - Controller structure ('P', 'PI', 'PD', 'PID')
    %   options   - Design options structure
    %   plantInfo - Plant information structure
    %
    % Outputs:
    %   K        - Emergency controller
    %   details  - Text description of the design process
    
    % Initialize details
    details = 'SPECIALIZED EMERGENCY CONTROLLER DESIGN\n';
    details = [details, '======================================\n'];
    details = [details, sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo))];
    
    % Extract key plant properties
    p = plantInfo.poles;
    z = plantInfo.zeros;
    
    % Extract unstable poles and analyze pattern
    unstable_poles = p(real(p) > 0);
    
    if length(unstable_poles) > 0
        details = [details, sprintf('\nPlant has %d unstable pole(s):\n', length(unstable_poles))];
        
        % Sort unstable poles by real part (most unstable first)
        [~, idx] = sort(real(unstable_poles), 'descend');
        unstable_poles = unstable_poles(idx);
        
        for i = 1:min(length(unstable_poles), 5)
            if imag(unstable_poles(i)) ~= 0
                details = [details, sprintf('  • %.4f + %.4fj\n', real(unstable_poles(i)), imag(unstable_poles(i)))];
            else
                details = [details, sprintf('  • %.4f\n', unstable_poles(i))];
            end
        end
    end
    
    % Extract RHP zeros and analyze pattern
    rhp_zeros = z(real(z) > 0);
    
    if length(rhp_zeros) > 0
        details = [details, sprintf('\nPlant has %d RHP zero(s):\n', length(rhp_zeros))];
        
        % Sort RHP zeros by real part (closest to imaginary axis first)
        [~, idx] = sort(real(rhp_zeros));
        rhp_zeros = rhp_zeros(idx);
        
        for i = 1:min(length(rhp_zeros), 5)
            if imag(rhp_zeros(i)) ~= 0
                details = [details, sprintf('  • %.4f + %.4fj\n', real(rhp_zeros(i)), imag(rhp_zeros(i)))];
            else
                details = [details, sprintf('  • %.4f\n', real(rhp_zeros(i)))];
            end
        end
    end
    
    % Define epsilon parameter for filtering
    if isfield(options, 'epsilon')
        epsilon = options.epsilon;
    else
        epsilon = 0.1;
    end
    
    % Increase filtering for robustness
    epsilon = max(epsilon, 0.2);
    
    % Use a multi-step approach for maximum robustness
    details = [details, '\n1. ATTEMPTING MULTI-PHASE EMERGENCY STABILIZATION:\n'];
    details = [details, '----------------------------------------------\n'];
    
    %% PHASE 1: Try direct pole cancellation with progressive gain reduction
    details = [details, 'PHASE 1: Direct unstable pole cancellation approach\n'];
    
    % First try direct pole cancellation approach
    try
        [K1, success1, details1] = directPoleCancellationApproach(G, unstable_poles, epsilon, rhp_zeros);
        details = [details, details1];
        
        if success1
            details = [details, 'Phase 1 successful! Direct cancellation approach stabilized the system.\n'];
            
            % Convert to requested structure if possible
            [K, success_convert] = convertToRequestedStructure(K1, G, structure, epsilon);
            
            if success_convert
                details = [details, sprintf('Successfully converted to requested %s structure.\n', structure)];
                return;
            else
                details = [details, 'Could not convert to requested structure while maintaining stability.\n'];
                details = [details, 'Using basic stabilizing controller instead.\n'];
                K = K1;
                return;
            end
        else
            details = [details, 'Phase 1 failed. Proceeding to Phase 2.\n\n'];
        end
    catch ME
        details = [details, sprintf('Error in Phase 1: %s\n', ME.message)];
        details = [details, 'Proceeding to Phase 2.\n\n'];
    end
    
    %% PHASE 2: Try layered approach (pre-stabilization + standard structure)
    details = [details, 'PHASE 2: Layered stabilization approach\n'];
    
    try
        [K2, success2, details2] = layeredStabilizationApproach(G, unstable_poles, structure, epsilon, rhp_zeros);
        details = [details, details2];
        
        if success2
            details = [details, 'Phase 2 successful! Layered approach stabilized the system.\n'];
            K = K2;
            return;
        else
            details = [details, 'Phase 2 failed. Proceeding to Phase 3.\n\n'];
        end
    catch ME
        details = [details, sprintf('Error in Phase 2: %s\n', ME.message)];
        details = [details, 'Proceeding to Phase 3.\n\n'];
    end
    
    %% PHASE 3: Try state-space approach if manageable order
    details = [details, 'PHASE 3: State-space stabilization approach\n'];
    
    try
        [K3, success3, details3] = stateSpaceStabilizationApproach(G, structure, epsilon, options, plantInfo);
        details = [details, details3];
        
        if success3
            details = [details, 'Phase 3 successful! State-space approach stabilized the system.\n'];
            K = K3;
            return;
        else
            details = [details, 'Phase 3 failed. Proceeding to Phase 4.\n\n'];
        end
    catch ME
        details = [details, sprintf('Error in Phase 3: %s\n', ME.message)];
        details = [details, 'Proceeding to Phase 4.\n\n'];
    end
    
    %% PHASE 4: Last resort - ultra-conservative approach
    details = [details, 'PHASE 4: Ultra-conservative stabilization approach\n'];
    
    try
        [K4, success4, details4] = ultraConservativeApproach(G, structure, epsilon, plantInfo);
        details = [details, details4];
        
        if success4
            details = [details, 'Phase 4 successful! Ultra-conservative approach stabilized the system.\n'];
            K = K4;
            return;
        else
            details = [details, 'Phase 4 failed. System may be extraordinarily difficult to control.\n'];
        end
    catch ME
        details = [details, sprintf('Error in Phase 4: %s\n', ME.message)];
    end
    
    %% Ultimate fallback - minimal controller
    details = [details, '\n2. ALL METHODS FAILED - CREATING MINIMAL CONTROLLER\n'];
    details = [details, '-----------------------------------------------\n'];
    details = [details, 'Creating minimal gain controller as absolute last resort\n'];
    
    % Create minimal controller based on structure
    switch structure
        case 'P'
            K = tf(0.00001, 1);
        case 'PI'
            K = tf([0.00001, 0.000001], [1, 0]);
        case 'PD'
            K = tf([0.0001, 0.00001], [0.1, 1]);
        case 'PID'
            K = tf([0.0001, 0.00001, 0.000001], [0.1, 1, 0]);
        otherwise
            K = tf(0.00001, 1);
    end
    
    details = [details, '\nWARNING: System appears to be extremely difficult to control.\n'];
    details = [details, 'Manual controller design is strongly recommended.\n'];
    details = [details, 'Consider redesigning the plant if possible.\n'];
end

%% Helper Functions

function [K, success, details] = directPoleCancellationApproach(G, unstable_poles, epsilon, rhp_zeros)
    % Direct pole cancellation approach
    
    details = '';
    success = false;
    
    % Check if there are unstable poles
    if isempty(unstable_poles)
        K = tf(0.1, 1);
        details = 'No unstable poles found. Using simple proportional controller.\n';
        
        % Check if stabilizes
        try
            T = feedback(G*K, 1);
            if all(real(pole(T)) < 0)
                success = true;
                details = [details, 'Simple proportional controller stabilizes the system.\n'];
                return;
            end
        catch
            % Do nothing, continue with more approaches
        end
        
        return;
    end
    
    % Initialize compensator
    num_K = 1;
    den_K = 1;
    
    % Create zeros to cancel each unstable pole
    details = [details, 'Creating cancellation zeros for unstable poles:\n'];
    
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
                
                details = [details, sprintf('  - Complex pole at %.4f+%.4fj\n', real_part, imag_part)];
                
                % Create zeros at unstable poles
                quad_term_num = [1, -2*real_part, real_part^2 + imag_part^2];
                
                % Create stable replacement poles far in LHP
                new_real_part = -abs(real_part) * 5;  % 5x damping
                new_quad_term = [1, -2*new_real_part, new_real_part^2 + imag_part^2];
                
                num_K = conv(num_K, quad_term_num);
                den_K = conv(den_K, new_quad_term);
            end
        else
            % For real poles
            details = [details, sprintf('  - Real pole at %.4f\n', pole_i)];
            
            % Create zero at unstable pole
            num_K = conv(num_K, [1, -pole_i]);
            
            % Create stable replacement pole far in LHP
            stable_pole_loc = -abs(pole_i) * 5;  % 5x margin for stability
            den_K = conv(den_K, [1, -stable_pole_loc]);
        end
    end
    
    % If RHP zeros exist, add high-frequency roll-off
    if ~isempty(rhp_zeros)
        details = [details, 'Plant has RHP zeros. Adding high-frequency roll-off.\n'];
        
        % Add roll-off based on closest RHP zero
        min_real_zero = min(real(rhp_zeros));
        rolloff_freq = min_real_zero * 0.2;  % Conservative limit
        
        rolloff_term = [1/rolloff_freq^2, sqrt(2)/rolloff_freq, 1];
        den_K = conv(den_K, rolloff_term);
        
        details = [details, sprintf('  - Added rolloff at %.4f rad/s\n', rolloff_freq)];
    end
    
    % Try progressively smaller gains until stability is achieved
    details = [details, 'Testing different gain levels for stability:\n'];
    
    % Initial gain set very conservatively
    base_gain = 0.001;
    
    for gain_scale = [1, 0.1, 0.01, 0.001, 0.0001, 0.00001]
        gain = base_gain * gain_scale;
        details = [details, sprintf('  - Testing gain = %.8f\n', gain)];
        
        % Create controller with current gain
        K_test = tf(gain * num_K, den_K);
        
        % Test closed-loop stability
        try
            T_test = feedback(G * K_test, 1);
            cl_poles = pole(T_test);
            
            if all(real(cl_poles) < 0)
                K = K_test;
                success = true;
                details = [details, sprintf('Success! Stable with gain = %.8f\n', gain)];
                
                % Show achieved stability margins if possible
                try
                    [Gm, Pm, Wcg, Wcp] = margin(G*K);
                    details = [details, sprintf('  - Phase margin: %.2f degrees\n', Pm)];
                    details = [details, sprintf('  - Gain margin: %.2f dB\n', 20*log10(Gm))];
                catch
                    % Skip margins if calculation fails
                end
                
                return;
            end
        catch ME
            details = [details, sprintf('  - Error testing gain: %s\n', ME.message)];
        end
    end
    
    % If we reach here, all gains failed
    K = tf(base_gain * 0.00001 * num_K, den_K);
    details = [details, 'All gain levels failed to stabilize the system.\n'];
end

function [K, success, details] = layeredStabilizationApproach(G, unstable_poles, structure, epsilon, rhp_zeros)
    % Layered approach: first stabilize, then add structure
    
    details = '';
    success = false;
    
    % Phase 2A: Create very basic stabilizing controller
    details = [details, 'Step 2A: Creating basic stabilizing controller\n'];
    
    % For each unstable pole region, create different stabilizing controller
    if isempty(unstable_poles)
        % No unstable poles - try simple high gain feedback
        K_stab = tf(1, 1);
        details = [details, 'No unstable poles. Using unity gain controller.\n'];
    elseif length(unstable_poles) == 1 && imag(unstable_poles(1)) == 0
        % Single real unstable pole - phase lead compensator
        p = real(unstable_poles(1));
        a = 0.1;  % Lead ratio
        T = 1/p;  % Time constant based on unstable pole
        
        num = [T, 1];
        den = [a*T, 1];
        
        K_stab = tf(0.1*num, den);  % Small initial gain
        
        details = [details, sprintf('Single real unstable pole. Using lead compensator with a=%.2f, T=%.4f\n', a, T)];
    elseif length(unstable_poles) > 0
        % Multiple or complex unstable poles - use pole placement approach
        % Create a custom loop-shaping controller
        
        % Start with the rightmost unstable pole
        p_max = unstable_poles(1);
        
        if imag(p_max) ~= 0
            % Complex unstable pole - use resonant compensator
            real_part = real(p_max);
            imag_part = imag(p_max);
            
            % Resonant numerator (creates notch at unstable frequency)
            num = [1, 0, imag_part^2];
            
            % Stable denominator (well-damped)
            den = [1, 2*abs(real_part), imag_part^2 + real_part^2];
            
            K_stab = tf(0.01*num, den);
            
            details = [details, sprintf('Complex unstable poles. Using resonant compensator.\n')];
        else
            % Multiple real unstable poles - use high-order lead compensator
            num = 1;
            den = 1;
            
            for i = 1:min(length(unstable_poles), 3)  % Handle up to 3 unstable poles
                p = real(unstable_poles(i));
                a = 0.1;  % Lead ratio
                T = 1/p;  % Time constant
                
                num_i = [T, 1];
                den_i = [a*T, 1];
                
                num = conv(num, num_i);
                den = conv(den, den_i);
            end
            
            K_stab = tf(0.001*num, den);
            
            details = [details, sprintf('Multiple real unstable poles. Using cascaded lead compensators.\n')];
        end
    end
    
    % Add rolloff for plants with RHP zeros
    if ~isempty(rhp_zeros)
        details = [details, 'Plant has RHP zeros. Adding rolloff to compensator.\n'];
        
        [num, den] = tfdata(K_stab, 'v');
        
        % Add rolloff based on closest RHP zero
        min_real_zero = min(real(rhp_zeros));
        rolloff_freq = min_real_zero * 0.2;  % Conservative
        
        rolloff_term = [1/rolloff_freq^2, sqrt(2)/rolloff_freq, 1];
        new_den = conv(den, rolloff_term);
        
        K_stab = tf(num, new_den);
    end
    
    % Try progressively smaller gains until stability is achieved
    details = [details, 'Testing base controller with different gains:\n'];
    
    [num, den] = tfdata(K_stab, 'v');
    stabilized = false;
    
    for gain_scale = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001]
        K_test = tf(num * gain_scale, den);
        details = [details, sprintf('  - Testing gain = %.6f\n', gain_scale)];
        
        try
            T_test = feedback(G * K_test, 1);
            cl_poles = pole(T_test);
            
            if all(real(cl_poles) < 0)
                K_stab = K_test;
                stabilized = true;
                details = [details, sprintf('Success! Base controller stabilizes with gain = %.6f\n', gain_scale)];
                
                % Show achieved stability margins if possible
                try
                    [Gm, Pm, Wcg, Wcp] = margin(G*K_stab);
                    details = [details, sprintf('  - Phase margin: %.2f degrees\n', Pm)];
                    details = [details, sprintf('  - Gain margin: %.2f dB\n', 20*log10(Gm))];
                catch
                    % Skip margins if calculation fails
                end
                
                break;
            end
        catch ME
            details = [details, sprintf('  - Error testing gain: %s\n', ME.message)];
        end
    end
    
    if ~stabilized
        details = [details, 'Base controller failed to stabilize system.\n'];
        K = K_stab;
        return;
    end
    
    % Phase 2B: Add requested structure on top of stabilizing controller
    details = [details, '\nStep 2B: Adding requested structure\n'];
    
    % Create stabilized plant
    G_stab = feedback(G * K_stab, 1);
    
    details = [details, sprintf('Created stabilized plant G_stab using base controller.\n')];
    
    % Now design standard controller for the stabilized plant
    try
        switch structure
            case 'P'
                % Simple P controller with very conservative gain
                K_struct = tf(0.1, 1);
                details = [details, 'Added P controller with gain = 0.1\n'];
                
            case 'PI'
                % PI with small integral gain to minimize disturbance
                Kp = 0.1;
                Ki = 0.01;
                K_struct = tf([Kp, Ki], [1, 0]);
                details = [details, sprintf('Added PI controller with Kp = %.3f, Ki = %.3f\n', Kp, Ki)];
                
            case 'PD'
                % PD with heavy filtering
                Kp = 0.1;
                Kd = 0.05;
                Tf = epsilon * 5;  % Extra filtering for robustness
                K_struct = tf([Kd, Kp], [Tf, 1]);
                details = [details, sprintf('Added PD controller with Kp = %.3f, Kd = %.3f, Tf = %.3f\n', Kp, Kd, Tf)];
                
            case 'PID'
                % PID with heavy filtering and small integral gain
                Kp = 0.1;
                Ki = 0.01;
                Kd = 0.05;
                Tf = epsilon * 5;  % Extra filtering for robustness
                K_struct = tf([Kd, Kp, Ki], [Tf, 1, 0]);
                details = [details, sprintf('Added PID controller with Kp = %.3f, Ki = %.3f, Kd = %.3f, Tf = %.3f\n', Kp, Ki, Kd, Tf)];
                
            otherwise
                % Default to P
                K_struct = tf(0.1, 1);
                details = [details, 'Added simple P controller with gain = 0.1\n'];
        end
        
        % Check if the structured controller stabilizes the stabilized plant
        try
            T_struct = feedback(G_stab * K_struct, 1);
            struct_poles = pole(T_struct);
            
            if all(real(struct_poles) < 0)
                details = [details, 'Structure controller stabilizes the stabilized plant!\n'];
                
                % Combine the controllers
                K_combined = series(K_struct, K_stab);
                
                % Final verification with original plant
                T_final = feedback(G * K_combined, 1);
                final_poles = pole(T_final);
                
                if all(real(final_poles) < 0)
                    K = K_combined;
                    success = true;
                    details = [details, 'Combined controller successfully stabilizes the original plant!\n'];
                    
                    % Show achieved stability margins if possible
                    try
                        [Gm, Pm, Wcg, Wcp] = margin(G*K);
                        details = [details, sprintf('Final controller phase margin: %.2f degrees\n', Pm)];
                        details = [details, sprintf('Final controller gain margin: %.2f dB\n', 20*log10(Gm))];
                    catch
                        % Skip margins if calculation fails
                    end
                    
                    return;
                else
                    details = [details, 'Combined controller fails to stabilize original plant.\n'];
                    details = [details, 'Using base stabilizing controller only.\n'];
                    K = K_stab;
                    success = true;
                    return;
                end
            else
                details = [details, 'Structure controller fails to stabilize the stabilized plant.\n'];
                details = [details, 'Using base stabilizing controller only.\n'];
                K = K_stab;
                success = true;
                return;
            end
        catch ME
            details = [details, sprintf('Error testing structured controller: %s\n', ME.message)];
            details = [details, 'Using base stabilizing controller only.\n'];
            K = K_stab;
            success = true;
            return;
        end
    catch ME
        details = [details, sprintf('Error creating structured controller: %s\n', ME.message)];
        details = [details, 'Using base stabilizing controller only.\n'];
        K = K_stab;
        success = true;
        return;
    end
end

function [K, success, details] = stateSpaceStabilizationApproach(G, structure, epsilon, options, plantInfo)
    % State-space based stabilization approach
    
    details = '';
    success = false;
    
    try
        % Get or create state-space model
        if isfield(options, 'stateSpace') && ~isempty(options.stateSpace)
            sys_ss = options.stateSpace;
            details = [details, 'Using provided state-space model.\n'];
        else
            % Try to create a state-space model
            try
                sys_ss = ss(G);
                details = [details, 'Successfully created state-space model from transfer function.\n'];
            catch ME1
                details = [details, sprintf('Basic state-space conversion failed: %s\n', ME1.message)];
                
                try
                    % Try more robust conversion methods
                    sys_ss = getEnhancedStateSpace(G, plantInfo);
                    details = [details, 'Successfully created enhanced state-space model.\n'];
                catch ME2
                    details = [details, sprintf('Enhanced conversion also failed: %s\n', ME2.message)];
                    details = [details, 'State-space approach cannot be used.\n'];
                    K = tf(0.001, 1);
                    return;
                end
            end
        end
        
        % Extract matrices
        A = sys_ss.A;
        B = sys_ss.B;
        C = sys_ss.C;
        D = sys_ss.D;
        
        % Verify controllability
        n = size(A, 1);
        
        if n > 20
            details = [details, 'System order too high for reliable state-space design.\n'];
            K = tf(0.001, 1);
            return;
        end
        
        Cm = ctrb(A, B);
        rank_C = rank(Cm);
        
        if rank_C < n
            details = [details, sprintf('System not fully controllable (rank %d < order %d).\n', rank_C, n)];
            
            if rank_C < n/2
                details = [details, 'Too many uncontrollable modes. State-space approach not suitable.\n'];
                K = tf(0.001, 1);
                return;
            else
                details = [details, 'Proceeding with reduced-order design.\n'];
            end
        else
            details = [details, 'System is controllable (good).\n'];
        end
        
        % Design state feedback controller
        details = [details, 'Designing state feedback controller with pole placement.\n'];
        
        % Create desired closed-loop poles
        p = eig(A);
        unstable_idx = real(p) >= 0;
        
        if all(~unstable_idx)
            details = [details, 'System appears stable in state-space form.\n'];
            details = [details, 'Using conservative pole placement for robustness.\n'];
            
            % Create conservative poles in LHP
            desired_poles = zeros(1, n);
            
            for i = 1:n
                desired_poles(i) = -1 - i;  % Simple stable pattern: -1, -2, -3, ...
            end
        else
            details = [details, 'Creating desired pole pattern for unstable system.\n'];
            
            % Create desired poles by reflecting unstable poles
            desired_poles = p;
            
            for i = 1:length(p)
                if real(p(i)) >= 0
                    if imag(p(i)) ~= 0
                        % For complex poles, maintain frequency but add damping
                        mag = abs(p(i));
                        desired_poles(i) = -mag * 0.5 + 1j * imag(p(i)) * 0.5;
                    else
                        % For real poles, reflect and move further left
                        desired_poles(i) = -abs(real(p(i))) * 2 - 1;
                    end
                else
                    % Move stable poles further left for robustness
                    desired_poles(i) = p(i) * 1.5;
                end
            end
            
            % Ensure complex poles come in conjugate pairs
            for i = 1:length(desired_poles)
                if imag(desired_poles(i)) ~= 0
                    conj_found = false;
                    
                    for j = 1:length(desired_poles)
                        if i ~= j && abs(desired_poles(i) - conj(desired_poles(j))) < 1e-6
                            conj_found = true;
                            break;
                        end
                    end
                    
                    if ~conj_found
                        % Make a conjugate pair
                        for j = 1:length(desired_poles)
                            if imag(desired_poles(j)) == 0
                                desired_poles(j) = conj(desired_poles(i));
                                break;
                            end
                        end
                    end
                end
            end
        end
        
        details = [details, 'Desired pole locations:\n'];
        for i = 1:min(length(desired_poles), 5)
            if imag(desired_poles(i)) ~= 0
                details = [details, sprintf('  • %.4f + %.4fj\n', real(desired_poles(i)), imag(desired_poles(i)))];
            else
                details = [details, sprintf('  • %.4f\n', desired_poles(i))];
            end
        end
        
        % Try multiple pole placement methods with error handling
        try
            K_state = place(A, B, desired_poles);
            details = [details, 'Successfully placed poles using place() function.\n'];
        catch ME1
            details = [details, sprintf('place() failed: %s\n', ME1.message)];
            
            try
                K_state = acker(A, B, desired_poles);
                details = [details, 'Successfully placed poles using acker() function.\n'];
            catch ME2
                details = [details, sprintf('acker() failed: %s\n', ME2.message)];
                
                try
                    % Last resort - simple design for dominant poles
                    K_state = zeros(1, n);
                    
                    % Use dot product with eigenvectors for basic control
                    [V, D] = eig(A);
                    lambda = diag(D);
                    
                    for i = 1:n
                        if real(lambda(i)) >= 0
                            vi = V(:, i);
                            bi = B' * vi;
                            
                            if abs(bi) > 1e-6
                                K_state = K_state - (real(lambda(i)) + 2) * vi' / bi;
                            end
                        end
                    end
                    
                    details = [details, 'Created basic feedback using eigenvector approach.\n'];
                catch ME3
                    details = [details, sprintf('All pole placement methods failed: %s\n', ME3.message)];
                    K = tf(0.001, 1);
                    return;
                end
            end
        end
        
        % Create state-space controller
        Ac = A - B*K_state;
        Bc = B;
        Cc = -K_state;
        Dc = 0;
        
        K_ss = ss(Ac, Bc, Cc, Dc);
        
        % Convert to transfer function
        K_tf = tf(K_ss);
        
        % Verify stabilization
        try
            T_ss = feedback(G * K_tf, 1);
            ss_poles = pole(T_ss);
            
            if all(real(ss_poles) < 0)
                details = [details, 'State-space controller successfully stabilizes the system!\n'];
                
                % Convert to requested structure if possible
                [K, success_convert] = convertToRequestedStructure(K_tf, G, structure, epsilon);
                
                if success_convert
                    details = [details, sprintf('Successfully converted to requested %s structure.\n', structure)];
                    success = true;
                    return;
                else
                    details = [details, 'Could not convert to requested structure while maintaining stability.\n'];
                    details = [details, 'Using state-space controller directly.\n'];
                    K = K_tf;
                    success = true;
                    return;
                end
            else
                details = [details, 'State-space controller fails to stabilize system.\n'];
                
                % Try scaling the gain
                [num, den] = tfdata(K_tf, 'v');
                
                for scale = [0.5, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001]
                    K_test = tf(num * scale, den);
                    
                    try
                        T_test = feedback(G * K_test, 1);
                        test_poles = pole(T_test);
                        
                        if all(real(test_poles) < 0)
                            details = [details, sprintf('Stabilized with gain scaling factor: %.4f\n', scale)];
                            K = K_test;
                            success = true;
                            return;
                        end
                    catch
                        % Skip this scale if analysis fails
                    end
                end
                
                details = [details, 'All scaling attempts failed.\n'];
                K = K_tf;
            end
        catch ME
            details = [details, sprintf('Error validating state-space controller: %s\n', ME.message)];
            K = K_tf;
        end
    catch ME
        details = [details, sprintf('Error in state-space approach: %s\n', ME.message)];
        K = tf(0.001, 1);
    end
end

function [K, success, details] = ultraConservativeApproach(G, structure, epsilon, plantInfo)
    % Ultra-conservative controller design as last resort
    
    details = '';
    success = false;
    
    % Get plant properties
    try
        p = plantInfo.poles;
        z = plantInfo.zeros;
    catch
        p = pole(G);
        try
            z = zero(G);
        catch
            z = [];
        end
    end
    
    % Create a series of extremely conservative controllers and test them
    controllers = {};
    
    % 1. Simple super-low gain
    controllers{1} = tf(0.0001, 1);
    
    % 2. Low-pass filter with extremely low bandwidth
    wc = 0.01;
    controllers{2} = tf(wc^2, [1, 2*wc, wc^2]);
    
    % 3. Controller with structure but tiny gains
    switch structure
        case 'P'
            controllers{3} = tf(0.0001, 1);
        case 'PI'
            controllers{3} = tf([0.0001, 0.00001], [1, 0]);
        case 'PD'
            controllers{3} = tf([0.0001, 0.0001], [0.1, 1]);
        case 'PID'
            controllers{3} = tf([0.0001, 0.0001, 0.00001], [0.1, 1, 0]);
        otherwise
            controllers{3} = tf(0.0001, 1);
    end
    
    % 4. Notch filter at RHP zero frequencies (if any)
    if ~isempty(z) && any(real(z) > 0)
        rhp_zeros = z(real(z) > 0);
        
        if length(rhp_zeros) > 0
            z_critical = rhp_zeros(1);
            if imag(z_critical) ~= 0
                w0 = abs(z_critical);
                zeta = 0.1;
                num = [1, 2*zeta*w0, w0^2];
                den = [1, 2*0.9*w0, w0^2];
                controllers{4} = tf(0.0001 * num, den);
            else
                w0 = real(z_critical) * 0.1;
                controllers{4} = tf([1, w0], [1, w0*10]);
            end
        end
    end
    
    % 5. Lead compensator with extreme parameters
    if ~isempty(p) && any(real(p) > 0)
        unstable_poles = p(real(p) > 0);
        if length(unstable_poles) > 0
            p_critical = unstable_poles(1);
            w0 = abs(p_critical) * 0.1;
            controllers{5} = tf([1, w0], [1, w0*0.01]);
        end
    end
    
    details = [details, 'Testing series of ultra-conservative controllers:\n'];
    
    % Try each controller
    for i = 1:length(controllers)
        if isempty(controllers{i})
            continue;
        end
        
        details = [details, sprintf('Testing controller type %d\n', i)];
        
        try
            T = feedback(G * controllers{i}, 1);
            cl_poles = pole(T);
            
            if all(real(cl_poles) < 0)
                details = [details, sprintf('Controller type %d successfully stabilizes the system!\n', i)];
                K = controllers{i};
                success = true;
                return;
            end
        catch ME
            details = [details, sprintf('Error testing controller type %d: %s\n', i, ME.message)];
        end
        
        % Try with scaled gains
        [num, den] = tfdata(controllers{i}, 'v');
        
        for scale = [0.1, 0.01, 0.001, 0.0001]
            K_test = tf(num * scale, den);
            
            try
                T_test = feedback(G * K_test, 1);
                test_poles = pole(T_test);
                
                if all(real(test_poles) < 0)
                    details = [details, sprintf('Controller type %d with gain scaling %.6f stabilizes the system!\n', i, scale)];
                    K = K_test;
                    success = true;
                    return;
                end
            catch
                % Skip this scale if analysis fails
            end
        end
    end
    
    % If all controllers fail, return the most basic one
    K = controllers{1};
    details = [details, 'All ultra-conservative controllers failed to stabilize the system.\n'];
end

function [K, success] = convertToRequestedStructure(K_original, G, structure, epsilon)
    % Convert a controller to the requested structure while maintaining stability
    
    success = false;
    
    % Get original controller response at key frequencies
    try
        w = logspace(-3, 3, 20);
        [mag_orig, phase_orig] = bode(K_original, w);
        mag_orig = squeeze(mag_orig);
        phase_orig = squeeze(phase_orig);
    catch
        % If frequency response fails, extract DC gain
        try
            dc_gain = dcgain(K_original);
            if isnan(dc_gain) || isinf(dc_gain)
                dc_gain = 1;
            end
        catch
            dc_gain = 1;
        end
    end
    
    % Create controller with requested structure that approximates original response
    try
        switch structure
            case 'P'
                % Simple P with gain matching at crossover or DC
                try
                    % Get gain at 1 rad/s or middle of frequency range
                    mid_idx = ceil(length(w)/2);
                    Kp = mag_orig(mid_idx);
                    
                    if isnan(Kp) || isinf(Kp) || Kp <= 0
                        Kp = 0.1;  % Default if calculation fails
                    end
                catch
                    Kp = 0.1;  % Default
                end
                
                K = tf(Kp, 1);
                
            case 'PI'
                % PI with gain and low-frequency response matching
                try
                    % Get phase and magnitude at key frequencies
                    low_idx = 1;                   % Lowest frequency
                    mid_idx = ceil(length(w)/2);   % Middle frequency
                    
                    % Set Kp to match gain at mid frequency
                    Kp = mag_orig(mid_idx) * 0.7;  % Reduced for stability
                    
                    % Set Ki to approximate integral action
                    if length(w) >= 3
                        % Estimate integral action from low-frequency phase
                        phase_diff = phase_orig(low_idx) - phase_orig(3);
                        if phase_diff > 30 && phase_diff < 150
                            % If phase suggests integral action
                            Ki = Kp * 0.1;
                        else
                            Ki = Kp * 0.01;  % Conservative
                        end
                    else
                        Ki = Kp * 0.01;  % Default conservative
                    end
                    
                    if isnan(Kp) || isinf(Kp) || Kp <= 0
                        Kp = 0.1;  % Default if calculation fails
                    end
                    
                    if isnan(Ki) || isinf(Ki) || Ki <= 0
                        Ki = 0.01;  % Default if calculation fails
                    end
                catch
                    Kp = 0.1;   % Default
                    Ki = 0.01;  % Default
                end
                
                K = tf([Kp, Ki], [1, 0]);
                
            case 'PD'
                % PD with gain and high-frequency response matching
                try
                    % Get phase and magnitude at key frequencies
                    mid_idx = ceil(length(w)/2);   % Middle frequency
                    high_idx = length(w);          % Highest frequency
                    
                    % Set Kp to match gain at mid frequency
                    Kp = mag_orig(mid_idx) * 0.7;  % Reduced for stability
                    
                    % Set Kd to approximate derivative action
                    if high_idx > mid_idx
                        % Estimate derivative action from high-frequency magnitude
                        mag_ratio = mag_orig(high_idx) / mag_orig(mid_idx);
                        if mag_ratio > 1.5
                            % If magnitude suggests derivative action
                            Kd = Kp * w(high_idx) / 10;
                        else
                            Kd = Kp * 0.1;  % Conservative
                        end
                    else
                        Kd = Kp * 0.1;  % Default conservative
                    end
                    
                    if isnan(Kp) || isinf(Kp) || Kp <= 0
                        Kp = 0.1;  % Default if calculation fails
                    end
                    
                    if isnan(Kd) || isinf(Kd) || Kd <= 0
                        Kd = 0.01;  % Default if calculation fails
                    end
                catch
                    Kp = 0.1;   % Default
                    Kd = 0.01;  % Default
                end
                
                Tf = epsilon;
                K = tf([Kd, Kp], [Tf, 1]);
                
            case 'PID'
                % PID with gain, low and high-frequency response matching
                try
                    % Get phase and magnitude at key frequencies
                    low_idx = 1;                   % Lowest frequency
                    mid_idx = ceil(length(w)/2);   % Middle frequency
                    high_idx = length(w);          % Highest frequency
                    
                    % Set Kp to match gain at mid frequency
                    Kp = mag_orig(mid_idx) * 0.7;  % Reduced for stability
                    
                    % Set Ki to approximate integral action
                    if length(w) >= 3
                        % Estimate integral action from low-frequency phase
                        phase_diff = phase_orig(low_idx) - phase_orig(3);
                        if phase_diff > 30 && phase_diff < 150
                            % If phase suggests integral action
                            Ki = Kp * 0.1;
                        else
                            Ki = Kp * 0.01;  % Conservative
                        end
                    else
                        Ki = Kp * 0.01;  % Default conservative
                    end
                    
                    % Set Kd to approximate derivative action
                    if high_idx > mid_idx
                        % Estimate derivative action from high-frequency magnitude
                        mag_ratio = mag_orig(high_idx) / mag_orig(mid_idx);
                        if mag_ratio > 1.5
                            % If magnitude suggests derivative action
                            Kd = Kp * w(high_idx) / 10;
                        else
                            Kd = Kp * 0.1;  % Conservative
                        end
                    else
                        Kd = Kp * 0.1;  % Default conservative
                    end
                    
                    if isnan(Kp) || isinf(Kp) || Kp <= 0
                        Kp = 0.1;  % Default if calculation fails
                    end
                    
                    if isnan(Ki) || isinf(Ki) || Ki <= 0
                        Ki = 0.01;  % Default if calculation fails
                    end
                    
                    if isnan(Kd) || isinf(Kd) || Kd <= 0
                        Kd = 0.01;  % Default if calculation fails
                    end
                catch
                    Kp = 0.1;   % Default
                    Ki = 0.01;  % Default
                    Kd = 0.01;  % Default
                end
                
                Tf = epsilon;
                K = tf([Kd, Kp, Ki], [Tf, 1, 0]);
                
            otherwise
                % Default to P controller if unknown structure
                K = tf(0.1, 1);
        end
        
        % Verify if the converted controller stabilizes the system
        T = feedback(G * K, 1);
        cl_poles = pole(T);
        
        if all(real(cl_poles) < 0)
            success = true;
        else
            % Try scaled gains for stability
            [num, den] = tfdata(K, 'v');
            
            for scale = [0.5, 0.2, 0.1, 0.05, 0.01]
                K_test = tf(num * scale, den);
                
                try
                    T_test = feedback(G * K_test, 1);
                    test_poles = pole(T_test);
                    
                    if all(real(test_poles) < 0)
                        K = K_test;
                        success = true;
                        return;
                    end
                catch
                    % Skip this scale if analysis fails
                end
            end
            
            % If all scales fail, return original
            K = K_original;
        end
    catch ME
        % If conversion fails, return original
        K = K_original;
    end
end

function sys_ss = getEnhancedStateSpace(G, plantInfo)
    % Create enhanced state-space representation with improved numerical properties
    
    % Try different approaches in sequence of increasing numerical robustness
    try
        % Try balanced realization first (typically good numerical properties)
        sys_ss = balreal(ss(G));
    catch
        try
            % Try model reduction if balreal fails
            order = length(pole(G));
            sys_ss = balred(ss(G), order);
        catch
            try
                % If model reduction fails, try canonical forms
                sys_ss = ss(G, 'canonical');
            catch
                % Last resort: manual construction from transfer function
                [num, den] = tfdata(G, 'v');
                
                % Check for ill-conditioning and apply regularization if needed
                if (max(abs(num)) / min(abs(num(num ~= 0))) > 1e10) || ...
                   (max(abs(den)) / min(abs(den(den ~= 0))) > 1e10)
                    % Scale coefficients to improve conditioning
                    scale = max(max(abs(num)), max(abs(den)));
                    num = num / scale;
                    den = den / scale;
                end
                
                % Handle potential issues with zeros at the end
                if abs(den(end)) < 1e-10
                    den = den(1:end-1);
                end
                
                if abs(num(end)) < 1e-10
                    num = num(1:end-1);
                end
                
                % Create manually using control canonical form
                n = length(den) - 1;  % System order
                
                % Control canonical form matrices
                A = zeros(n);
                A(1:n-1, 2:n) = eye(n-1);
                A(n, :) = -den(2:end) ./ den(1);
                
                B = zeros(n, 1);
                B(n) = 1 ./ den(1);
                
                C = zeros(1, n);
                if length(num) <= n
                    C(1:length(num)) = num ./ den(1);
                else
                    C = num(2:n+1) ./ den(1);
                end
                
                D = 0;
                if length(num) > n
                    D = num(1) ./ den(1);
                end
                
                sys_ss = ss(A, B, C, D);
            end
        end
    end
end

function info_str = getPlantInfoString(plantInfo)
    % Create a formatted string with plant information
    
    % Initialize output string
    info_str = '';
    
    % Add stability information
    if plantInfo.isUnstable
        info_str = [info_str, 'Unstable, '];
    else
        info_str = [info_str, 'Stable, '];
    end
    
    % Add integrator information
    if plantInfo.hasIntegrator
        info_str = [info_str, 'Has integrator, '];
    end
    
    % Add RHP zeros information
    if plantInfo.hasRHPZeros
        info_str = [info_str, 'Non-minimum phase, '];
    end
    
    % Add delay information
    if plantInfo.hasDelay
        info_str = [info_str, 'Has delay, '];
    end
    
    % Add order information
    if plantInfo.isHighOrder
        info_str = [info_str, 'High-order, '];
    else
        info_str = [info_str, 'Low-order, '];
    end
    
    % Add DC gain if available
    if ~isnan(plantInfo.dcGain) && ~isinf(plantInfo.dcGain)
        info_str = [info_str, sprintf('DC gain=%.3g', plantInfo.dcGain)];
    else
        info_str = [info_str, 'Infinite DC gain'];
    end
end