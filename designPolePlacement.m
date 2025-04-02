function [K, details] = designPolePlacement(G, structure, options, plantInfo)
% DESIGNPOLEPLACEMENT Controller design using pole placement method
%
% Implements the pole placement (Zustandsrückführung mit Polvorgabe) method
% as described in section 6.6 of the Regelungstechnik document. Enhanced to work
% as an automated controller design method with more robust error handling.
%
% Inputs:
%   G        - Plant transfer function
%   structure - Controller structure ('P', 'PI', 'PD', 'PID')
%   options   - Structure with design parameters:
%     .bandwidth - Desired bandwidth in rad/s (default: 1)
%     .damping   - Desired damping ratio (default: 0.8)
%     .epsilon   - Filter parameter for D-term (default: 0.1)
%     .stateSpace - Optional state-space model (if provided from outside)
%     .stateMatrices - Optional {A, B, C, D} cell array
%   plantInfo - Structure with plant analysis information
%
% Outputs:
%   K       - Designed controller as transfer function
%   details - Text description with design details

    % Start with detailed information about the method
    details = 'Pole Placement Design Method (Zustandsrückführung mit Polvorgabe)\n';
    details = [details, '--------------------------------------------------------\n'];
    details = [details, sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo))];

    % Default values for bandwidth and damping
    if ~isfield(options, 'bandwidth')
        options.bandwidth = 1;
    end
    
    if ~isfield(options, 'damping')
        options.damping = 0.8;
    end
    
    if ~isfield(options, 'epsilon')
        options.epsilon = 0.1;
    end
    
    % Extract bandwidth and damping
    omega = options.bandwidth;
    zeta = options.damping;
    epsilon = options.epsilon;
    
    details = [details, sprintf('Desired bandwidth: %.4f rad/s\n', omega)];
    details = [details, sprintf('Desired damping ratio: %.4f\n', zeta)];
    details = [details, sprintf('Derivative filter coefficient: %.4f\n', epsilon)];
    
    % First, get or convert to state space
    try
        % Check if state-space model is provided
        if isfield(options, 'stateSpace') && ~isempty(options.stateSpace)
            sys_ss = options.stateSpace;
            details = [details, 'Using provided state-space model.\n'];
        elseif isfield(options, 'stateMatrices') && ~isempty(options.stateMatrices)
            A = options.stateMatrices{1};
            B = options.stateMatrices{2};
            C = options.stateMatrices{3};
            D = options.stateMatrices{4};
            sys_ss = ss(A, B, C, D);
            details = [details, 'Using provided state-space matrices.\n'];
        else
            % Try the direct conversion first
            sys_ss = ss(G);
            details = [details, 'Successfully converted to state-space representation.\n'];
        end
    catch ME
        % If direct conversion fails, try a different approach for problematic cases
        details = [details, sprintf('Direct state-space conversion failed: %s\n', ME.message)];
        details = [details, 'Attempting alternative conversion approach...\n'];
        
        try
            % First check if there's a pole at the origin (integrator)
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
                sys_ss = ss(G_reduced);
            else
                % Try to balance the system to improve numerical conditioning
                [num, den] = tfdata(G, 'v');
                sys_tf = tf(num, den);
                sys_ss = balreal(ss(sys_tf));
            end
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
                    % Do nothing, keep k = 1
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
                sys_ss = ss(A, B, C, D);
                
                details = [details, 'Using numerically balanced model.\n'];
            catch ME3
                details = [details, sprintf('All conversion attempts failed: %s\n', ME3.message)];
                error('Cannot convert to workable state-space model for pole placement.');
            end
        end
    end
    
    % Calculate desired pole locations based on bandwidth and damping
    if plantInfo.isUnstable
        % For unstable systems, we need more conservative pole placement
        % Scale down bandwidth to avoid over-aggressive control
        effective_omega = omega * 0.3;
        
        details = [details, sprintf('Unstable system detected. Reducing effective bandwidth to %.4f rad/s\n', effective_omega)];
        
        % For unstable systems, prefer real poles (more robust)
        p1 = -effective_omega;
        p2 = -effective_omega * 2;
        desired_poles = [p1, p2];
        
        % Add more poles for higher-order systems
        plant_order = length(pole(G));
        if plant_order > 2
            extra_poles = -effective_omega * (3:plant_order+1);
            desired_poles = [desired_poles, extra_poles];
        end
    else
        % For stable systems, use standard second-order pole placement
        if zeta < 1
            % Underdamped - complex conjugate poles
            real_part = -zeta * omega;
            imag_part = omega * sqrt(1 - zeta^2);
            p1 = complex(real_part, imag_part);
            p2 = complex(real_part, -imag_part);
            desired_poles = [p1, p2];
        else
            % Critically damped or overdamped - real poles
            p1 = -omega;
            p2 = -omega * (2*zeta - 1);  % Second pole further left for overdamped
            desired_poles = [p1, p2];
        end
        
        % For higher-order systems, add more poles further left
        plant_order = length(pole(G));
        if plant_order > 2
            extra_poles = -omega * (3:plant_order+1);
            desired_poles = [desired_poles, extra_poles];
        end
    end
    
    details = [details, 'Desired pole locations:\n'];
    for i = 1:length(desired_poles)
        if imag(desired_poles(i)) ~= 0
            details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(desired_poles(i)), imag(desired_poles(i)))];
        else
            details = [details, sprintf('  p%d = %.4f\n', i, real(desired_poles(i)))];
        end
    end
    
    % Extract state-space matrices
    A = sys_ss.A;
    B = sys_ss.B;
    C = sys_ss.C;
    D = sys_ss.D;
    
    % Check controllability
    Co = ctrb(A, B);
    controllability_rank = rank(Co);
    system_order = size(A, 1);
    
    details = [details, sprintf('Controllability check: rank %d out of %d\n', controllability_rank, system_order)];
    
    if controllability_rank < system_order
        details = [details, 'WARNING: System is not fully controllable!\n'];
        details = [details, 'Using partial pole placement for controllable subspace.\n'];
        
        % For non-controllable systems, we would need more advanced techniques
        % Here we'll use a simpler approach with regularization
        A = A + 1e-6 * eye(size(A));
        Co = ctrb(A, B);
        controllability_rank = rank(Co);
        details = [details, sprintf('After regularization: rank %d out of %d\n', controllability_rank, system_order)];
    end
    
    try
        % Attempt pole placement using place() function
        K_state = place(A, B, desired_poles);
        details = [details, 'Successfully computed state feedback gain K.\n'];
        details = [details, sprintf('K = [%s]\n', mat2str(K_state, 4))];
        
    catch ME
        details = [details, sprintf('Standard pole placement failed: %s\n', ME.message)];
        details = [details, 'Using robust fallback approach...\n'];
        
        try
            % Instead of precise pole placement, use a robust pole region approach
            % This creates a more conservative controller that's more likely to work
            
            % Get plant information
            A = sys_ss.A;
            B = sys_ss.B;
            C = sys_ss.C;
            n = size(A, 1);
            
            % For unstable plants, focus on stabilization first
            if plantInfo.isUnstable
                % Find unstable poles
                p = eig(A);
                unstable_idx = find(real(p) > 0);
                
                % Create a stabilizing gain by placing unstable poles in LHP
                p_desired = p;
                for i = 1:length(unstable_idx)
                    p_desired(unstable_idx(i)) = -abs(real(p(unstable_idx(i))));
                end
                
                % Try to compute a gain that just stabilizes the system
                K_state = zeros(1, n);
                
                % Simple pole shifting approach for each unstable pole
                for i = 1:length(unstable_idx)
                    shift_vector = zeros(n, 1);
                    shift_vector(unstable_idx(i)) = 1;
                    
                    % Compute how much to shift this pole
                    shift_amount = real(p(unstable_idx(i))) * 2;
                    
                    % Update the gain
                    K_state = K_state + shift_amount * shift_vector' * pinv(B);
                end
                
                details = [details, 'Created stabilizing gain using pole shifting.\n'];
            else
                % For stable plants, use a simple output feedback approach
                K_state = -options.bandwidth * pinv(B) * ones(n, 1);
                details = [details, 'Created simple feedback gain.\n'];
            end
            
            details = [details, sprintf('K = [%s]\n', mat2str(K_state, 4))];
        catch ME2
            details = [details, sprintf('Fallback approach failed: %s\n', ME2.message)];
            
            % Last resort: use a very simple gain
            n = size(A, 1);
            K_state = ones(1, n);
            
            % Try to scale the gain for better results
            try
                K_state = K_state * sign(B(1));
            catch
                % If that fails, just use the unscaled version
            end
            
            details = [details, 'Using simple unity gain as last resort.\n'];
        end
    end
    
    % Now convert the state feedback controller to the requested structure
    details = [details, sprintf('\nConverting state feedback to %s controller...\n', structure)];
    
    % Calculate the closed-loop system with state feedback
    Acl = A - B * K_state;
    sys_cl = ss(Acl, B, C, D);
    
    % Create controller based on the required structure
    switch structure
        case 'P'
            % Extract a proportional gain that approximates the state feedback
            Kp = abs(K_state(1));  % Use absolute value for stability
            if Kp < 0.1
                Kp = 0.1; % Minimum gain for stability
            elseif Kp > 100
                Kp = 100; % Maximum gain to prevent excessive control effort
            end
            
            K = tf(Kp, 1);
            details = [details, sprintf('Extracted P controller: Kp = %.4f\n', Kp)];
            
        case 'PI'
            % For PI, we need to include an integrator
            Kp = abs(K_state(1));
            
            % Relate Ki to bandwidth and Kp
            if plantInfo.isUnstable
                Ki = Kp * omega / 20;  % Very small Ki for unstable systems
            else
                Ki = Kp * omega / 5;   % Larger Ki for stable systems
            end
            
            % Limit values to prevent instability
            if Kp < 0.1
                Kp = 0.1;
            elseif Kp > 100
                Kp = 100;
            end
            
            if Ki < 0.01
                Ki = 0.01;
            elseif Ki > 50
                Ki = 50;
            end
            
            K = tf([Kp, Ki], [1, 0]);
            details = [details, sprintf('Designed PI controller: Kp = %.4f, Ki = %.4f\n', Kp, Ki)];
            
        case 'PD'
            % For PD, extract the derivative action
            if length(K_state) >= 2
                Kp = abs(K_state(1));
                Kd = abs(K_state(2)) / omega;  % Scale by bandwidth
            else
                Kp = abs(K_state(1));
                Kd = Kp / omega;  % Relate Kd to bandwidth
            end
            
            % Limit values to prevent instability
            if Kp < 0.1
                Kp = 0.1;
            elseif Kp > 100
                Kp = 100;
            end
            
            if Kd < 0.01
                Kd = 0.01;
            elseif Kd > 50
                Kd = 50;
            end
            
            % Add filtering for derivative term
            K = tf([Kd, Kp], [epsilon*Kd, 1]);
            details = [details, sprintf('Designed PD controller: Kp = %.4f, Kd = %.4f\n', Kp, Kd)];
            
        case 'PID'
            % For PID, extract all three actions
            Kp = abs(K_state(1));
            
            % Computing Ki and Kd based on the state feedback and desired dynamics
            if plantInfo.isUnstable
                Ki = Kp * omega / 20;  % Very small Ki for unstable systems
            else
                Ki = Kp * omega / 5;   % Larger Ki for stable systems
            end
            
            if length(K_state) >= 2
                Kd = abs(K_state(2)) / omega;  % Scale by bandwidth
            else
                Kd = Kp / omega^2;  % Relate Kd to bandwidth
            end
            
            % Limit values to prevent instability
            if Kp < 0.1
                Kp = 0.1;
            elseif Kp > 100
                Kp = 100;
            end
            
            if Ki < 0.01
                Ki = 0.01;
            elseif Ki > 50
                Ki = 50;
            end
            
            if Kd < 0.01
                Kd = 0.01;
            elseif Kd > 50
                Kd = 50;
            end
            
            % Add filtering for derivative term
            K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
            details = [details, sprintf('Designed PID controller: Kp = %.4f, Ki = %.4f, Kd = %.4f\n', Kp, Ki, Kd)];
            
        otherwise
            error('Unsupported controller structure for pole placement method');
    end
    
    % Verify the closed-loop system with the designed controller
    try
        closed_loop = feedback(G*K, 1);
        cl_poles = pole(closed_loop);
        
        % Check stability of the result
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
        
        % If unstable, adjust the controller
        if ~is_stable && ~strcmpi(structure, 'P')
            details = [details, '\nAttempting to stabilize controller by adjusting parameters...\n'];
            
            % Reduce integral gain and increase derivative gain
            switch structure
                case 'PI'
                    [num, den] = tfdata(K, 'v');
                    Kp = num(1);
                    Ki = num(2);
                    Ki_new = Ki * 0.1;  % Significantly reduce integral gain
                    K = tf([Kp, Ki_new], [1, 0]);
                    details = [details, sprintf('Adjusted PI controller: Kp = %.4f, Ki = %.4f\n', Kp, Ki_new)];
                    
                case 'PD'
                    [num, den] = tfdata(K, 'v');
                    Kd = num(1);
                    Kp = num(2);
                    Kd_new = Kd * 3;  % Significantly increase derivative gain
                    K = tf([Kd_new, Kp], [epsilon*Kd_new, 1]);
                    details = [details, sprintf('Adjusted PD controller: Kp = %.4f, Kd = %.4f\n', Kp, Kd_new)];
                    
                case 'PID'
                    [num, den] = tfdata(K, 'v');
                    Kd = num(1);
                    Kp = num(2);
                    Ki = num(3);
                    Ki_new = Ki * 0.05;  % Drastically reduce integral gain
                    Kd_new = Kd * 5;   % Significantly increase derivative gain
                    K = tf([Kd_new, Kp, Ki_new], [epsilon*Kd_new, 1, 0]);
                    details = [details, sprintf('Adjusted PID controller: Kp = %.4f, Ki = %.4f, Kd = %.4f\n', Kp, Ki_new, Kd_new)];
            end
            
            % Check if adjustment succeeded
            try
                closed_loop_adj = feedback(G*K, 1);
                cl_poles_adj = pole(closed_loop_adj);
                is_stable_adj = all(real(cl_poles_adj) < 0);
                
                if is_stable_adj
                    details = [details, 'Success! Adjusted controller stabilizes the system.\n'];
                    
                    details = [details, 'New closed-loop poles:\n'];
                    for i = 1:length(cl_poles_adj)
                        if imag(cl_poles_adj(i)) ~= 0
                            details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(cl_poles_adj(i)), imag(cl_poles_adj(i)))];
                        else
                            details = [details, sprintf('  p%d = %.4f\n', i, real(cl_poles_adj(i)))];
                        end
                    end
                else
                    details = [details, 'Warning: System remains unstable after adjustment.\n'];
                    details = [details, 'Consider using a more robust design method for this plant.\n'];
                    
                    % Last resort: make a very conservative controller
                    if strcmpi(structure, 'PID')
                        Kp_last = 0.1;
                        Ki_last = 0.001;
                        Kd_last = 10;
                        K = tf([Kd_last, Kp_last, Ki_last], [epsilon*Kd_last, 1, 0]);
                        details = [details, sprintf('Last resort PID: Kp = %.4f, Ki = %.4f, Kd = %.4f\n', Kp_last, Ki_last, Kd_last)];
                        
                        % Final check
                        try
                            closed_loop_last = feedback(G*K, 1);
                            cl_poles_last = pole(closed_loop_last);
                            is_stable_last = all(real(cl_poles_last) < 0);
                            
                            if is_stable_last
                                details = [details, 'Last resort controller is stable!\n'];
                            else
                                details = [details, 'Even last resort controller is unstable. This system is extremely challenging.\n'];
                            end
                        catch
                            details = [details, 'Could not verify last resort controller stability.\n'];
                        end
                    end
                end
            catch
                details = [details, 'Could not verify adjusted controller stability.\n'];
            end
        end
        
        % Add a performance evaluation
        try
            % Evaluate step response if system is stable
            if is_stable
                [y, t] = step(closed_loop, linspace(0, 10/omega, 1000));
                stepInfo = stepinfo(y, t);
                
                details = [details, '\nClosed-loop performance:\n'];
                
                if isfield(stepInfo, 'RiseTime') && ~isnan(stepInfo.RiseTime)
                    details = [details, sprintf('Rise time: %.4f s\n', stepInfo.RiseTime)];
                end
                
                if isfield(stepInfo, 'SettlingTime') && ~isnan(stepInfo.SettlingTime)
                    details = [details, sprintf('Settling time: %.4f s\n', stepInfo.SettlingTime)];
                end
                
                if isfield(stepInfo, 'Overshoot') && ~isnan(stepInfo.Overshoot)
                    details = [details, sprintf('Overshoot: %.2f%%\n', stepInfo.Overshoot)];
                end
                
                % Add stability margins
                try
                    [Gm, Pm, ~, ~] = margin(G*K);
                    Gm_dB = 20*log10(Gm);
                    
                    details = [details, sprintf('Gain margin: %.2f dB\n', Gm_dB)];
                    details = [details, sprintf('Phase margin: %.2f degrees\n', Pm)];
                    
                    % Evaluate controller quality based on margins
                    score = 0;
                    
                    % Base score from stability
                    if is_stable
                        score = 60;  % Starting score for stable system
                    else
                        score = 0;   % Unstable system gets 0
                    end
                    
                    % Add points for good margins
                    if Gm_dB > 10
                        score = score + 10;
                    elseif Gm_dB > 6
                        score = score + 5;
                    end
                    
                    if Pm > 60
                        score = score + 10;
                    elseif Pm > 45
                        score = score + 5;
                    end
                    
                    % Add points for good time response
                    if isfield(stepInfo, 'Overshoot') && ~isnan(stepInfo.Overshoot)
                        if stepInfo.Overshoot < 5
                            score = score + 10;
                        elseif stepInfo.Overshoot < 15
                            score = score + 5;
                        elseif stepInfo.Overshoot > 25
                            score = score - 5;
                        end
                    end
                    
                    % Adjust for complexity of the plant
                    if plantInfo.isUnstable || plantInfo.hasRHPZeros
                        score = score + 5;  % Bonus for handling difficult plants
                    end
                    
                    details = [details, sprintf('\nController Score: %.2f/100\n', score)];
                catch
                    details = [details, 'Could not calculate stability margins.\n'];
                end
            end
        catch
            details = [details, 'Could not analyze closed-loop performance.\n'];
        end
    catch ME
        details = [details, sprintf('\nError in closed-loop analysis: %s\n', ME.message)];
    end
    
    return;
end
