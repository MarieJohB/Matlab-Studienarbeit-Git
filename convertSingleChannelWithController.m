function [G_tf, K_tf] = convertSingleChannelWithController(ssModel, K_full, inputIdx, outputIdx)
    % CONVERTSINGLECONTROLCHANNEL Creates transfer functions for a single input-output channel
    % with its corresponding controller
    %
    % Parameters:
    %   ssModel - State space model (ss object)
    %   K_full - Full state feedback gain matrix
    %   inputIdx - Selected input index
    %   outputIdx - Selected output index
    %
    % Returns:
    %   G_tf - Plant transfer function for the selected channel
    %   K_tf - Corresponding controller transfer function

    % Extract state-space matrices
    A = ssModel.A;
    B = ssModel.B;
    C = ssModel.C;
    D = ssModel.D;
    Ts = ssModel.Ts;
    
    % Get system dimensions
    [n_states, n_inputs] = size(B);
    [n_outputs, ~] = size(C);
    
    % Validate input and output indices
    if inputIdx > n_inputs || inputIdx < 1
        error('Input index out of range');
    end
    if outputIdx > n_outputs || outputIdx < 1
        error('Output index out of range');
    end
    
    % Check K matrix dimensions
    [k_rows, k_cols] = size(K_full);
    if k_cols ~= n_states || k_rows ~= n_inputs
        warning('K matrix dimensions (%dx%d) don''t match expected dimensions (%dx%d)', ...
            k_rows, k_cols, n_inputs, n_states);
    end
    
    % 1. Create plant transfer function for the selected channel
    try
        % Convert from state space to transfer function for specific I/O pair
        G_tf = tf(ssModel, inputIdx, outputIdx);
        
        % Simplify and clean up numerical noise
        [num_g, den_g] = tfdata(G_tf, 'v');
        num_g = cleanup_tf_coefficients(num_g);
        den_g = cleanup_tf_coefficients(den_g);
        G_tf = tf(num_g, den_g, Ts);
        
        disp(['Created plant transfer function G for input ' num2str(inputIdx) ...
             ' to output ' num2str(outputIdx)]);
    catch ME
        error('Error creating plant transfer function: %s', ME.message);
    end
    
    % 2. Create a corresponding controller transfer function
    try
        % Extract the relevant row from the K matrix for the selected input
        K_row = K_full(inputIdx, :);
        
        % Method 1: Direct Observer-Based Approach
        % For this input-output pair, create a dynamic controller by augmenting
        % the system to include an observer for estimating states
        % We'll use a Luenberger observer with poles 4x faster than closed-loop poles
        
        % Calculate closed-loop poles
        cl_poles = eig(A - B(:,inputIdx) * K_row);
        
        % Calculate observer poles (4x faster)
        obs_poles = cl_poles * 4;
        
        % Design observer gain
        try
            L = place(A', C(outputIdx,:)', obs_poles)';
        catch
            % If place fails, try using a simpler approach
            L = 4 * ones(n_states, 1);
            warning('Observer design with place failed, using simplified observer');
        end
        
        % Create observer-based controller
        % This implements K(s) = K * (sI - A + BK + LC)^-1 * L
        Ac = A - B(:,inputIdx)*K_row - L*C(outputIdx,:);
        Bc = L;
        Cc = K_row;
        Dc = 0;
        
        % Create transfer function of the controller
        K_tf = ss(Ac, Bc, Cc, Dc);
        K_tf = tf(K_tf);
        
        % Clean up numerical noise in the controller
        [num_k, den_k] = tfdata(K_tf, 'v');
        num_k = cleanup_tf_coefficients(num_k);
        den_k = cleanup_tf_coefficients(den_k);
        K_tf = tf(num_k, den_k, Ts);
        
        disp(['Created observer-based controller for input ' num2str(inputIdx) ...
             ' and state feedback gain K']);
    catch ME
        warning('Error creating dynamic controller: %s', ME.message);
        % Fallback: Use simpler proportional-only controller
        K_tf = tf(norm(K_row), 1);
        warning('Using simplified controller with gain K = %.4f', norm(K_row));
    end
end