% Helper function to design a stabilizing controller for unstable plants
function K_tf = designStabilizingController(G_tf, Ts)
    % This function designs a controller specifically focused on stabilizing an unstable plant
    
    % Extract plant info
    [num_g, den_g] = tfdata(G_tf, 'v');
    
    % Get plant order
    plant_order = length(den_g) - 1;
    
    % Try different controller structures based on instability characteristics
    
    % First try: PID-type controller
    % For unstable plants, PID controllers often need to be carefully tuned
    if Ts > 0  % Discrete-time
        % Discrete PID parameters
        kp = 1.0;
        ki = 0.2;
        kd = 0.2;
        
        % Simple discrete approximation of PID
        pid_num = [kp+ki+kd, -(kp+2*kd), kd];
        pid_den = [1, -1, 0];
        
        K_tf = tf(pid_num, pid_den, Ts);
    else  % Continuous-time
        % PID parameters
        kp = 1.0;
        ki = 0.1;
        kd = 0.1;
        
        % Create PID with proper filtering
        pid_num = [kd, kp, ki];
        pid_den = [0.01, 1, 0];  % Small time constant for derivative term
        
        K_tf = tf(pid_num, pid_den);
    end
    
    % Check if this stabilizes the system
    T_tf = feedback(series(G_tf, K_tf), 1);
    
    if Ts > 0  % Discrete-time
        is_stable = all(abs(pole(T_tf)) < 1);
    else  % Continuous-time
        is_stable = all(real(pole(T_tf)) < 0);
    end
    
    % If PID worked, return it
    if is_stable
        disp('PID-type controller successfully stabilizes the system');
        return;
    end
    
    % Second try: Lead compensator for phase improvement
    if Ts > 0  % Discrete-time
        % Discrete lead compensator
        a = 0.1;  % a < 1 for lead
        T = 1.0;  % Time constant
        lead_num = [1, -a];
        lead_den = [1, -1/T];
        
        K_tf = tf(lead_num, lead_den, Ts);
    else  % Continuous-time
        % Continuous lead compensator
        a = 0.1;  % a < 1 for lead
        T = 1.0;  % Time constant
        lead_num = [T, 1];
        lead_den = [a*T, 1];
        
        K_tf = tf(lead_num, lead_den);
    end
    
    % Check stability
    T_tf = feedback(series(G_tf, K_tf), 1);
    
    if Ts > 0  % Discrete-time
        is_stable = all(abs(pole(T_tf)) < 1);
    else  % Continuous-time
        is_stable = all(real(pole(T_tf)) < 0);
    end
    
    % If lead compensator worked, return it
    if is_stable
        disp('Lead compensator successfully stabilizes the system');
        return;
    end
    
    % Third try: More aggressive notch + lead design to target unstable poles
    p = pole(G_tf);
    
    % Find unstable poles
    if Ts > 0  % Discrete-time
        unstable_poles = p(abs(p) >= 1);
    else  % Continuous-time
        unstable_poles = p(real(p) >= 0);
    end
    
    % If no unstable poles found but system is unstable, use a default approach
    if isempty(unstable_poles)
        % Default fallback - high-gain stabilization
        disp('Using high-gain approach for stabilization');
        K_tf = tf(10, 1, Ts);
        return;
    end
    
    % Design a compensator to specifically target the dominant unstable pole
    [~, idx] = max(real(unstable_poles));
    dominant_pole = unstable_poles(idx);
    
    if Ts > 0  % Discrete-time
        % Create compensator that adds a zero near the unstable pole
        % and a pole inside the unit circle
        zero_loc = dominant_pole;
        pole_loc = 0.5 * dominant_pole / abs(dominant_pole);  % Inside unit circle
        
        K_num = [1, -zero_loc];
        K_den = [1, -pole_loc];
        
        % Add high-frequency roll-off
        K_num = conv(K_num, [1, 0]);
        K_den = conv(K_den, [1, 0.5]);
        
        K_tf = tf(K_num, K_den, Ts);
    else  % Continuous-time
        % Create compensator that adds a zero near the unstable pole
        % and a stable pole
        zero_loc = dominant_pole;
        pole_loc = -2 * abs(real(dominant_pole));  % Stable pole
        
        if imag(zero_loc) ~= 0
            % For complex poles, create complex conjugate pair
            K_num = conv([1, -zero_loc], [1, -conj(zero_loc)]);
            K_den = conv([1, -pole_loc], [1, -2*pole_loc]);
        else
            % For real poles
            K_num = [1, -zero_loc];
            K_den = [1, -pole_loc];
        end
        
        % Clean up to real coefficients if there was numerical error
        K_num = real(K_num);
        K_den = real(K_den);
        
        K_tf = tf(K_num, K_den);
    end
    
    % Add gain to improve response
    K_tf = 5 * K_tf;
    
    % Check one more time
    T_tf = feedback(series(G_tf, K_tf), 1);
    
    if Ts > 0  % Discrete-time
        is_stable = all(abs(pole(T_tf)) < 1);
    else  % Continuous-time
        is_stable = all(real(pole(T_tf)) < 0);
    end
    
    if is_stable
        disp('Custom pole-targeting compensator successfully stabilizes the system');
    else
        disp('Failed to stabilize with specialized compensator, using basic gain controller');
        K_tf = tf(5, 1, Ts);  % Last resort: simple gain
    end
end