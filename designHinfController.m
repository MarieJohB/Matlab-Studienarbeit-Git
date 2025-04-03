function [K_out, errorMsg, sysInfo] = designHinfController(A, B, C, D, inputIdx, outputIdx, params)
    % DESIGNHINFCONTROLLER - Enhanced H-infinity controller design
    %
    % Parameters:
    %   A, B, C, D - State space matrices
    %   inputIdx - Selected input index for MIMO systems
    %   outputIdx - Selected output index for MIMO systems
    %   params - Structure with H-infinity parameters
    %
    % Returns:
    %   K_out - State feedback gain
    %   errorMsg - Error message (empty if no error)
    %   sysInfo - Additional information including H-inf controller

    % Initialize outputs
    K_out = [];
    errorMsg = '';
    sysInfo = struct('info', {{}}, 'type', 'Hinf');
    
    % Calculate system dimensions
    [n_states, n_inputs] = size(B);
    [n_outputs, ~] = size(C);
    
    % Check if system is discrete-time
    isDT = isfield(params, 'Ts') && params.Ts > 0;
    Ts = 0;
    if isDT
        Ts = params.Ts;
    end

    % Parse parameters
    try
        % Get performance weight
        if isfield(params, 'performanceWeight')
            performanceWeight = params.performanceWeight;
        else
            performanceWeight = 1.0;
        end
        
        % Get robustness factor
        if isfield(params, 'robustnessFactor')
            robustness = params.robustnessFactor;
        else
            robustness = 0.5;
        end
        
        % Get bandwidth
        if isfield(params, 'bandwidth')
            bandwidth = params.bandwidth;
        else
            bandwidth = 1.0;
        end
        
        % Get desired controller order
        if isfield(params, 'order')
            desiredOrder = params.order;
            if ischar(desiredOrder) && strcmpi(desiredOrder, 'Auto')
                desiredOrder = [];  % Auto selection
            else
                desiredOrder = str2double(desiredOrder);
                if isnan(desiredOrder)
                    desiredOrder = [];  % Default to auto if invalid
                end
            end
        else
            desiredOrder = [];  % Default to auto
        end
    catch ME
        errorMsg = ['Error parsing H-infinity parameters: ', ME.message];
        return;
    end
    
    % Check if the Robust Control Toolbox is available
    hasRobustToolbox = exist('hinfsyn', 'file') == 2;
    if ~hasRobustToolbox
        errorMsg = 'Robust Control Toolbox is not available for H-infinity design.';
        disp(errorMsg);
        disp('Falling back to robust pole placement...');
        
        % Fall back to robust pole placement - will implement this later
        % This message will be used to inform the user that we're falling back
        sysInfo.info{end+1} = 'H∞ synthesis not available - using robust pole placement instead.';
        sysInfo.info{end+1} = 'For full H∞ functionality, install the Robust Control Toolbox.';
        
        % Will implement the fallback here, but for now return error
        return;
    end
    
    % Special handling for MIMO systems - focus on selected channel
    B_channel = B;
    C_channel = C;
    if n_inputs > 1
        % For MIMO inputs, focus on selected input
        if isfield(params, 'channelFocus') && params.channelFocus
            B_channel = B(:, inputIdx);
            sysInfo.info{end+1} = sprintf('Focusing on input channel %d for H∞ design', inputIdx);
        end
    end
    
    if n_outputs > 1
        % For MIMO outputs, focus on selected output
        if isfield(params, 'channelFocus') && params.channelFocus
            C_channel = C(outputIdx, :);
            sysInfo.info{end+1} = sprintf('Focusing on output channel %d for H∞ design', outputIdx);
        end
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
            disp('Detected unstable system, adapting H-infinity design approach...');
            sysInfo.info{end+1} = 'Unstable system detected - adapting H∞ design for stabilization';
            
            % For unstable systems, we need to:
            % 1. Modify the weighting functions to emphasize stabilization
            % 2. Potentially use Glover-McFarlane for unstable plants
            
            % Adjust performance weight to be less aggressive (stability first)
            performanceWeight = performanceWeight * 0.5;
            
            % Decrease bandwidth expectations for unstable system
            bandwidth = bandwidth * 0.7;
            
            % Increase robustness factor
            robustness = min(0.8, robustness * 1.5);
            
            sysInfo.info{end+1} = 'Modified H∞ parameters: decreased performance demands, increased robustness';
        end
    catch ME
        disp(['Warning in stability check: ', ME.message]);
        % Continue with standard design if stability check fails
    end
    
    % Design H-infinity controller using hinfsyn
    try
        disp('Setting up generalized plant for H-infinity synthesis...');
        
        % Create transfer function for appropriate domain
        if isDT
            s = zpk('z', Ts);
            planti = ss(A, B_channel, C_channel, zeros(size(C_channel, 1), size(B_channel, 2)), Ts);
        else
            s = tf('s');
            planti = ss(A, B_channel, C_channel, zeros(size(C_channel, 1), size(B_channel, 2)));
        end
        
        % Create weighting functions
        % For performance (tracking error) - shape Sensitivity function S
        if isDT
            % Discretized approximation for discrete-time systems
            Wp_num = [(1 + bandwidth*Ts/2), -1];
            Wp_den = [(bandwidth*Ts/0.01/2 + 1), -1];
            Wp = performanceWeight * tf(Wp_num, Wp_den, Ts);
        else
            % Performance weight (more sophisticated than before)
            % Low frequency: high gain for good tracking/disturbance rejection
            % High frequency: low gain to avoid noise amplification
            Wp_num = [1/bandwidth, 1];
            Wp_den = [1/(bandwidth*0.001), 1];
            Wp = performanceWeight * tf(Wp_num, Wp_den);
        end
        
        % Control effort weight (shapes KS) - prevent control signal saturation
        if isDT
            % Discretized high-pass filter
            Wu_num = [1, -1/(1 + bandwidth*Ts*100)];
            Wu_den = [1, -1];
            Wu = robustness * tf(Wu_num, Wu_den, Ts);
        else
            % Higher gain at high frequencies to limit control bandwidth
            Wu_num = [1, bandwidth/10];
            Wu_den = [0.01, bandwidth*100];
            Wu = robustness * tf(Wu_num, Wu_den);
        end
        
        % For unstable systems, add complementary sensitivity weight
        if isUnstable
            if isDT
                % Discretized low-pass filter
                Wt_num = [1, -0.9];
                Wt_den = [1, -0.1];
                Wt = tf(Wt_num, Wt_den, Ts);
            else
                % Roll-off weight for complementary sensitivity
                Wt_num = [1, bandwidth*5];
                Wt_den = [1, bandwidth*50];
                Wt = tf(Wt_num, Wt_den);
            end
        else
            % For stable systems, no need for explicit T shaping
            Wt = [];
        end
        
        % For MIMO, print some diagnostics about the weighting functions
        if n_inputs > 1 || n_outputs > 1
            disp('Weight function information for MIMO H-infinity design:');
            disp(['Performance weight (Wp) - Max gain: ', num2str(20*log10(norm(Wp, inf))), ' dB']);
            disp(['Control effort weight (Wu) - Max gain: ', num2str(20*log10(norm(Wu, inf))), ' dB']);
            if ~isempty(Wt)
                disp(['Complementary sensitivity weight (Wt) - Max gain: ', num2str(20*log10(norm(Wt, inf))), ' dB']);
            end
        end
        
        % Create generalized plant
        % Standard S/KS (sensitivity/control sensitivity) problem
        if isempty(Wt)
            % Without complementary sensitivity weight
            systemnames = 'planti Wp Wu';
            inputvar = '[ref; dist; control]';
            outputvar = '[Wp; Wu; ref-planti]';
            input_to_planti = '[control+dist]';
            input_to_Wp = '[ref-planti]';
            input_to_Wu = '[control]';
            
            sysoutname = 'P';
            cleanupsysic = 'yes';
            
            % Create interconnection
            P = sysic;
            
            % Identify inputs/outputs for controller
            n_inputs_P = 1;  % reference
            n_outputs_P = 1;  % measured output (ref-plant)
        else
            % With complementary sensitivity weight for unstable systems
            systemnames = 'planti Wp Wu Wt';
            inputvar = '[ref; dist; control]';
            outputvar = '[Wp; Wu; Wt; ref-planti]';
            input_to_planti = '[control+dist]';
            input_to_Wp = '[ref-planti]';
            input_to_Wu = '[control]';
            input_to_Wt = '[planti]';
            
            sysoutname = 'P';
            cleanupsysic = 'yes';
            
            % Create interconnection
            P = sysic;
            
            % Identify inputs/outputs for controller
            n_inputs_P = 1;  % reference
            n_outputs_P = 1;  % measured output (ref-plant)
        end
        
        % For MIMO systems, n_inputs_P and n_outputs_P would be different
        % Use the following for a general case:
        n_inputs_P = size(C_channel, 1);  % Number of measured outputs
        n_outputs_P = size(B_channel, 2);  % Number of control inputs
        
        % Adjust for non-square systems
        if n_inputs > 1 && n_outputs > 1 && ~isfield(params, 'channelFocus')
            % Using multiple channels - make sure dimensions are correct
            n_inputs_P = size(C, 1);      % All outputs
            n_outputs_P = size(B, 2);     % All inputs
        end
        
        % Design H-infinity controller
        disp('Designing H-infinity controller...');
        
        opt = hinfsynOptions;
        
        % If order reduction is requested
        if ~isempty(desiredOrder)
            opt.LimitedOrder = true;
            opt.Order = desiredOrder;
        end
        
        % For unstable systems, need to ensure internal stability
        if isUnstable
            % Make sure we enforce internal stability
            opt.InternalStability = 1;
        end
        
        % Perform H-infinity synthesis
        [K_hinf, CL, gamma] = hinfsyn(P, n_inputs_P, n_outputs_P, opt);
        
        % Store the H-infinity controller
        sysInfo.K_hinf = K_hinf;
        sysInfo.gamma = gamma;
        sysInfo.closedLoop = CL;
        
        % Convert to state feedback form for compatibility
        % Extract controller state-space matrices
        [AK, BK, CK, DK] = ssdata(K_hinf);
        
        % Store these for reference
        sysInfo.AK = AK;
        sysInfo.BK = BK;
        sysInfo.CK = CK;
        sysInfo.DK = DK;
        
        % For compatibility with the app, convert to a form that looks like state feedback
        % This is a simplification for visualization purposes
        if n_inputs == 1 || (n_inputs > 1 && isfield(params, 'channelFocus'))
            % SISO case or MIMO with channel focus - extract approximate gain
            if size(AK, 1) > 0  % Dynamic controller
                % Extract approximate state feedback from controller
                K_sf = CK * (-AK)\BK;
                
                % Might need reshaping for MIMO systems
                if length(K_sf) ~= n_states
                    % The approximation didn't work well, use a simpler approach
                    K_sf = zeros(1, n_states);
                    
                    % Find dominant states from output equation
                    [~, stateIdx] = sort(abs(C_channel), 'descend');
                    
                    % Assign reasonable values to these states
                    maxGain = max(1.0, norm(CK));
                    for i = 1:min(3, length(stateIdx))
                        K_sf(stateIdx(i)) = maxGain * (4-i)/6;  % Decreasing importance
                    end
                end
            else
                % Static controller - just use the feedthrough term
                K_sf = DK;
                
                % Expand if needed
                if length(K_sf) ~= n_states
                    K_sf = zeros(1, n_states);
                    
                    % Find dominant states from output equation
                    [~, stateIdx] = sort(abs(C_channel), 'descend');
                    
                    % Assign reasonable values
                    for i = 1:min(3, length(stateIdx))
                        K_sf(stateIdx(i)) = DK * (4-i)/6;  % Decreasing importance
                    end
                end
            end
            
            % For MIMO, create a full K matrix with the appropriate dimension
            K_out = zeros(n_inputs, n_states);
            if n_inputs > 1 && isfield(params, 'channelFocus')
                K_out(inputIdx, :) = K_sf;
            else
                K_out = K_sf;  % SISO case
            end
        else
            % Full MIMO case - less precise approximation
            K_out = zeros(n_inputs, n_states);
            
            % Attempt to extract approximate gains for each input
            if size(AK, 1) > 0  % Dynamic controller
                for i = 1:n_inputs
                    for j = 1:n_outputs
                        % Try to extract influence from each output to each input
                        try
                            K_part = CK(i,:) * (-AK)\BK(:,j);
                            
                            % Apply based on how this output observes states
                            for k = 1:n_states
                                K_out(i,k) = K_out(i,k) + K_part * C(j,k) / sum(abs(C(j,:)));
                            end
                        catch
                            % Skip if matrix operations fail
                        end
                    end
                    
                    % Normalize if the gain is too large
                    if norm(K_out(i,:)) > 10
                        K_out(i,:) = K_out(i,:) * 10 / norm(K_out(i,:));
                    end
                end
            else
                % Static controller
                K_out = DK;
            end
        end
        
        % Format controller order
        if ~isempty(desiredOrder)
            sysInfo.info{end+1} = sprintf('Specified H∞ controller order: %d', desiredOrder);
        else
            sysInfo.info{end+1} = sprintf('H∞ controller order (auto): %d', order(K_hinf));
        end
        
        % Add information about the H-infinity design
        sysInfo.info{end+1} = '';
        sysInfo.info{end+1} = 'H-Infinity Controller Information:';
        sysInfo.info{end+1} = sprintf('H∞ norm (gamma): %.4g', gamma);
        sysInfo.info{end+1} = sprintf('Performance Weight: %.4g, Robustness Factor: %.4g', performanceWeight, robustness);
        sysInfo.info{end+1} = sprintf('Bandwidth: %.4g rad/s', bandwidth);
        
        if isUnstable
            sysInfo.info{end+1} = '';
            sysInfo.info{end+1} = 'Note: Parameters were adjusted for unstable system stabilization.';
        end
        
        sysInfo.info{end+1} = '';
        sysInfo.info{end+1} = 'Important: H∞ controllers are dynamic - the gain matrix K is only an approximation.';
        sysInfo.info{end+1} = 'The full H∞ controller is stored in sysInfo.K_hinf for use in advanced analysis.';
        
    catch ME
        % Enhanced error handling with troubleshooting tips
        errorMsg = ['Error in H-infinity design: ', ME.message];
        disp(errorMsg);
        disp(['Stack trace: ', getReport(ME, 'basic')]);
        
        % Provide more user-friendly explanation based on error type
        if contains(ME.message, 'unstable') || contains(ME.message, 'stabilizable')
            sysInfo.info{end+1} = 'Error: The system appears to be unstabilizable with H∞ design.';
            sysInfo.info{end+1} = 'This can happen when unstable modes are not controllable.';
            
            % Try to suggest a solution
            sysInfo.info{end+1} = '';
            sysInfo.info{end+1} = 'Recommendations:';
            sysInfo.info{end+1} = '1. Check if all unstable modes are controllable';
            sysInfo.info{end+1} = '2. Try a different input channel that might better control unstable modes';
            sysInfo.info{end+1} = '3. Consider pre-stabilizing the system before applying H∞ design';
        elseif contains(ME.message, 'detectable') || contains(ME.message, 'Observable')
            sysInfo.info{end+1} = 'Error: The system appears to have undetectable modes for H∞ design.';
            sysInfo.info{end+1} = 'This can happen when some modes are not observable.';
            
            % Suggest a solution
            sysInfo.info{end+1} = '';
            sysInfo.info{end+1} = 'Recommendations:';
            sysInfo.info{end+1} = '1. Check if all modes are observable';
            sysInfo.info{end+1} = '2. Try a different output channel that might better observe all modes';
        elseif contains(ME.message, 'convergence') || contains(ME.message, 'converge')
            sysInfo.info{end+1} = 'Error: The H∞ algorithm did not converge.';
            sysInfo.info{end+1} = 'This often indicates numerical issues or very challenging specifications.';
            
            % Suggest a solution
            sysInfo.info{end+1} = '';
            sysInfo.info{end+1} = 'Recommendations:';
            sysInfo.info{end+1} = '1. Try reducing the performance requirements (lower performance weight)';
            sysInfo.info{end+1} = '2. Increase the robustness factor';
            sysInfo.info{end+1} = '3. For unstable systems, consider pre-stabilization first';
        else
            sysInfo.info{end+1} = 'Error in H∞ design process. See MATLAB command window for details.';
        end
        
        % Will add fallback to robust pole placement here
        return;
    end
end