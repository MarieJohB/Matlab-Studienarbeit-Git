% Helper function to create an observer-based controller
function K_tf = createObserverBasedController(A, B, C, K, inputIdx, outputIdx, Ts)
    % Extract relevant matrices for the channel
    B_channel = B(:, inputIdx);
    C_channel = C(outputIdx, :);
    
    % Extract the relevant row of K
    if size(K, 1) >= inputIdx
        K_row = K(inputIdx, :);
    else
        K_row = K(1, :);  % Use first row as fallback
    end
    
    % Create closed-loop poles for observer design
    cl_poles = eig(A - B_channel * K_row);
    
    % Calculate observer poles (typically faster than closed-loop)
    obs_poles = cl_poles * 3;  % Observer 3x faster
    
    % Limit pole speed to avoid numerical issues
    max_speed = 100;
    for i = 1:length(obs_poles)
        if Ts == 0  % Continuous-time
            if real(obs_poles(i)) > 0  % If unstable, flip to stable
                obs_poles(i) = -abs(real(obs_poles(i))) + 1j*imag(obs_poles(i));
            elseif abs(real(obs_poles(i))) > max_speed  % If too fast, limit speed
                angle_i = angle(obs_poles(i));
                obs_poles(i) = -max_speed * exp(1j*angle_i);
            end
        else  % Discrete-time
            if abs(obs_poles(i)) > 1  % If unstable, scale to unit circle
                obs_poles(i) = obs_poles(i) / (1.1 * abs(obs_poles(i)));
            end
        end
    end
    
    % Design observer gain using pole placement
    try
        L = place(A', C_channel', obs_poles)';
    catch ME
        disp(['Observer pole placement failed: ' ME.message]);
        disp('Using alternative observer design approach');
        
        % Alternative approach using DARE/CARE
        if Ts > 0  % Discrete-time
            [~, L] = dlqe(A, eye(size(A,1)), C_channel, eye(size(A,1)), eye(size(C_channel,1)));
        else  % Continuous-time
            [~, L] = lqe(A, eye(size(A,1)), C_channel, eye(size(A,1)), eye(size(C_channel,1)));
        end
    end
    
    % Create observer-based controller
    Ac = A - B_channel*K_row - L*C_channel;
    Bc = L;
    Cc = K_row;
    Dc = 0;
    
    % Create state-space model and convert to transfer function
    obs_ss = ss(Ac, Bc, Cc, Dc, Ts);
    K_tf = tf(obs_ss);
    
    disp('Created observer-based controller');
end