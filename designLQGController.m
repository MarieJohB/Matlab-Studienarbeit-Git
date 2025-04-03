function [K_out, errorMsg, sysInfo] = designLQGController(A, B, C, D, inputIdx, outputIdx, params)
    % DESIGNLQGCONTROLLER - Enhanced LQG controller design with full observer dynamics
    %
    % Parameters:
    %   A, B, C, D - State space matrices
    %   inputIdx - Selected input index for MIMO systems
    %   outputIdx - Selected output index for MIMO systems
    %   params - Structure with LQG parameters
    %
    % Returns:
    %   K_out - State feedback gain
    %   errorMsg - Error message (empty if no error)
    %   sysInfo - Additional information including observer gain

    % Initialize outputs
    K_out = [];
    errorMsg = '';
    sysInfo = struct('info', {{}}, 'type', 'LQG');
    
    % Calculate system dimensions
    [n_states, n_inputs] = size(B);
    [n_outputs, ~] = size(C);
    
    % Check if system is discrete-time
    isDT = isfield(params, 'Ts') && params.Ts > 0;

    % Parse Q and R matrices from parameters
    try
        if isfield(params, 'qDiag') && ~isempty(params.qDiag)
            qDiag = params.qDiag;
            Q = diag(qDiag);
        else
            % Default Q
            Q = eye(n_states);
        end
        
        if isfield(params, 'rDiag') && ~isempty(params.rDiag)
            rDiag = params.rDiag;
            R = diag(rDiag);
        else
            % Default R
            R = eye(n_inputs);
        end
        
        % Process and measurement noise covariances for Kalman filter
        if isfield(params, 'qnValue') && ~isempty(params.qnValue)
            qnValue = params.qnValue;
            Qn = eye(n_states) * qnValue;
        else
            % Default process noise
            Qn = eye(n_states) * 0.1;
        end
        
        if isfield(params, 'rnValue') && ~isempty(params.rnValue)
            rnValue = params.rnValue;
            Rn = eye(n_outputs) * rnValue;
        else
            % Default measurement noise
            Rn = eye(n_outputs) * 0.1;
        end
    catch ME
        errorMsg = ['Error parsing LQG parameters: ', ME.message];
        return;
    end
    
    % Special handling for unstable systems
    try
        % Check stability
        if isDT
            isUnstable = any(abs(eig(A)) >= 1);
        else
            isUnstable = any(real(eig(A)) >= 0);
        end
        
        if isUnstable
            disp('Detected unstable system, modifying LQG design approach...');
            
            % For unstable systems, do modal analysis to identify and target unstable modes
            [V, D] = eig(A);
            eigenvalues = diag(D);
            
            % Find unstable modes
            if isDT
                unstableIdx = find(abs(eigenvalues) >= 1);
            else
                unstableIdx = find(real(eigenvalues) >= 0);
            end
            
            % Increase penalty for states involved in unstable modes
            for i = 1:length(unstableIdx)
                % Get the eigenvector for this unstable mode
                modeIdx = unstableIdx(i);
                v = V(:, modeIdx);
                
                % Find states with significant contribution to this mode
                [~, maxStates] = sort(abs(v), 'descend');
                
                % Increase Q values for states that contribute most to unstable modes
                for j = 1:min(3, length(maxStates))  % Focus on top 3 contributing states
                    stateIdx = maxStates(j);
                    Q(stateIdx, stateIdx) = Q(stateIdx, stateIdx) * 10;  % Higher penalty
                end
            end
            
            % Make Kalman filter more aggressive for unstable systems
            Qn = Qn * 5;  % Higher process noise makes filter more responsive
            Rn = Rn * 0.5;  % Lower measurement noise makes filter trust measurements more
            
            disp('Adjusted Q, Qn, and Rn for unstable system.');
            sysInfo.info{end+1} = 'Unstable system detected - modified LQG parameters for stabilization';
        end
    catch ME
        disp(['Warning in stability check: ', ME.message]);
        % Continue with standard design if stability check fails
    end
    
    % Select appropriate B and C matrices for MIMO systems
    if n_inputs > 1 || n_outputs > 1
        disp('MIMO system detected in LQG design');
        
        % If we're focusing on a specific channel but using full-system LQG:
        if isfield(params, 'channelFocus') && params.channelFocus
            sysInfo.info{end+1} = sprintf('Focusing LQG design on input %d to output %d channel', inputIdx, outputIdx);
            
            % For MIMO LQG with channel focus, adjust R to emphasize selected input
            if n_inputs > 1
                for i = 1:n_inputs
                    if i ~= inputIdx
                        R(i,i) = R(i,i) * 10;  % Penalize other inputs more
                    end
                end
            end
            
            % Similarly adjust Rn to emphasize selected output
            if n_outputs > 1
                for i = 1:n_outputs
                    if i ~= outputIdx
                        Rn(i,i) = Rn(i,i) * 2;  % Trust other outputs less
                    end
                end
            end
        end
    end
    
    % Design LQG controller - first try with standard functions
    try
        disp('Designing LQR component...');
        
        % Design the LQR controller
        if isDT
            [K_lqr, S, e] = dlqr(A, B, Q, R);
        else
            [K_lqr, S, e] = lqr(A, B, Q, R);
        end
        
        disp('Designing Kalman filter component...');
        
        % Design the Kalman filter
        try
            % Try using the kalman function if available
            kalmansys = ss(A, [eye(n_states), B], C, [zeros(n_outputs, n_states), D]);
            [kest, L, P] = kalman(kalmansys, Qn, Rn);
            sysInfo.info{end+1} = 'Used MATLAB kalman function for filter design';
        catch ME
            disp(['Kalman function error: ', ME.message]);
            disp('Falling back to DARE approach for Kalman filter design');
            
            % Fall back to solving the filter Riccati equation directly
            if isDT
                [P, ~, ~] = dare(A', C', Qn, Rn);
            else
                [P, ~, ~] = care(A', C', Qn, Rn);
            end
            
            L = P * C' / (C * P * C' + Rn);
            
            if isDT
                L = A * L;  % Adjust for discrete-time case
            end
            
            sysInfo.info{end+1} = 'Used Riccati equation solving for Kalman filter design';
        end
        
        % Store the designed gains
        K_out = K_lqr;
        sysInfo.K_lqr = K_lqr;
        sysInfo.L = L;
        sysInfo.P = P;
        
        % Create full LQG controller model
        sysInfo.Ac = [A-B*K_lqr, B*K_lqr; zeros(n_states, n_states), A-L*C];
        sysInfo.Bc = [zeros(n_states, n_outputs); L];
        sysInfo.Cc = [K_lqr, zeros(n_inputs, n_states)];
        sysInfo.Dc = zeros(n_inputs, n_outputs);
        
        % Create a state-space model of the full LQG controller
        if isDT
            sysInfo.K_ss = ss(A-L*C-B*K_lqr, L, K_lqr, zeros(n_inputs, n_outputs), params.Ts);
        else
            sysInfo.K_ss = ss(A-L*C-B*K_lqr, L, K_lqr, zeros(n_inputs, n_outputs));
        end
                
        % Store the closed-loop eigenvalues
        sysInfo.closedLoopEigenvalues = e;
        sysInfo.info{end+1} = sprintf('LQR closed-loop eigenvalues: %s', formatPoles(e));
        
        % Add information about the controller design
        sysInfo.info{end+1} = '';
        sysInfo.info{end+1} = 'LQG Controller Information:';
        sysInfo.info{end+1} = sprintf('Q diagonal: [%s]', formattedArrayStr(diag(Q)));
        sysInfo.info{end+1} = sprintf('R diagonal: [%s]', formattedArrayStr(diag(R)));
        sysInfo.info{end+1} = sprintf('Process noise (Qn): %.4g, Measurement noise (Rn): %.4g', qnValue, rnValue);
        sysInfo.info{end+1} = '';
        sysInfo.info{end+1} = 'Important: LQG is a dynamic controller - the gain matrix K is only part of the controller.';
        sysInfo.info{end+1} = 'The full controller includes a Kalman filter/observer that estimates states before feedback.';
        
    catch ME
        % Enhanced error handling with context
        errorMsg = ['Error in LQG design: ', ME.message];
        disp(errorMsg);
        
        % Try to return a simplified controller if possible
        try
            % Design simple LQR without Kalman filter
            if isDT
                [K_lqr, ~, ~] = dlqr(A, B, Q, R);
            else
                [K_lqr, ~, ~] = lqr(A, B, Q, R);
            end
            
            K_out = K_lqr;
            sysInfo.info{end+1} = '';
            sysInfo.info{end+1} = 'WARNING: Full LQG design failed, using LQR only.';
            sysInfo.info{end+1} = sprintf('Error details: %s', ME.message);
            sysInfo.info{end+1} = 'Observer (Kalman filter) could not be designed.';
        catch ME2
            errorMsg = ['Error in fallback LQR design: ', ME2.message];
            return;
        end
    end
end

% Helper function to format poles for display
function str = formatPoles(poles)
    % Start with opening bracket
    str = '[';
    
    % Process each pole
    for i = 1:length(poles)
        if i > 1
            str = [str, ', '];
        end
        
        % Format based on whether the pole is real or complex
        if imag(poles(i)) == 0
            str = [str, sprintf('%.4g', real(poles(i)))];
        else
            str = [str, sprintf('%.4g %+.4gi', real(poles(i)), imag(poles(i)))];
        end
    end
    
    % Close the bracket
    str = [str, ']'];
end
