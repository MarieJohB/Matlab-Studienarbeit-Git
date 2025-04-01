function [K, details] = design_controller_auto(G, method, structure, options)
% DESIGN_CONTROLLER_AUTO Automatic controller design with various methods
% Enhanced version that prioritizes state-space representations when available
% 
% Inputs:
%   G         - Plant model (transfer function or state-space)
%   method    - Design method (string)
%   structure - Controller structure: 'P', 'PI', 'PD', 'PID'
%   options   - Structure with optional parameters:
%     .epsilon     - Filter parameter for D-term (default: 0.1)
%     .phaseMargin - Desired phase margin in degrees (default: 45)
%     .bandwidth   - Desired bandwidth in rad/s (default: 1)
%     .settlingTime - Desired settling time in s (default: 5)
%     .robustness  - Robustness level: 'Low', 'Medium', 'High' (default: 'Medium')
%     .overshoot   - Desired overshoot in % (default: 10)
%     .goal        - Optimization goal: 'Tracking', 'Disturbance Rejection', 'Robustness' (default: 'Tracking')
%
% Outputs:
%   K        - Designed controller as transfer function
%   details  - Text description with design details

    % Default values for missing options
    if nargin < 4
        options = struct();
    end
    
    if ~isfield(options, 'epsilon')
        options.epsilon = 0.1;
    end
    
    if ~isfield(options, 'phaseMargin')
        options.phaseMargin = 45;
    end
    
    if ~isfield(options, 'bandwidth')
        options.bandwidth = 1;
    end
    
    if ~isfield(options, 'settlingTime')
        options.settlingTime = 5;
    end
    
    if ~isfield(options, 'robustness')
        options.robustness = 'Medium';
    end
    
    if ~isfield(options, 'overshoot')
        options.overshoot = 10;
    end
    
    if ~isfield(options, 'goal')
        options.goal = 'Tracking';
    end
    
    if ~isfield(options, 'damping')
        options.damping = 0.8;
    end
    
    % Determine if input is state-space or transfer function
    isStateSpace = isa(G, 'ss');
    
    % Store raw state-space model if available for methods that can use it directly
    if isStateSpace
        options.stateSpace = G;  % Store original state-space model
        [A, B, C, D] = ssdata(G);
        options.stateMatrices = {A, B, C, D}; % Store matrices separately
    end
    
    % Pre-analyze the plant to determine characteristics
    plantInfo = analyzePlant(G);
    
    % For highly unstable systems, adjust parameters
    if plantInfo.isUnstable
        % Get a measure of instability
        unstable_poles = plantInfo.poles(real(plantInfo.poles) > 0);
        max_real_part = max(real(unstable_poles));
        
        % For highly unstable systems, use more conservative settings
        if max_real_part > 5 || length(unstable_poles) > 1
            disp('System is highly unstable - adjusting control parameters for better stability');
            
            % Reduce bandwidth for very unstable systems 
            if ~isfield(options, 'userSetBandwidth') || ~options.userSetBandwidth
                options.bandwidth = min(options.bandwidth, max_real_part/10);
                disp(['Adjusted bandwidth to ', num2str(options.bandwidth), ' rad/s']);
            end
            
            % Increase damping for highly unstable systems
            options.damping = min(1.5, options.damping * 1.5);
            
            % Set high robustness
            options.robustness = 'High';
        end
    end
    
    % Select design method
    try
        switch method
            case 'Ziegler-Nichols (Oscillation)'
                [K, details] = designZieglerNicholsOscillation(G, structure, options.epsilon, plantInfo);
                
            case 'Ziegler-Nichols (Step)'
                [K, details] = designZieglerNicholsStep(G, structure, options.epsilon, plantInfo);
                
            case 'Aström'
                [K, details] = designAstrom(G, structure, options.epsilon, plantInfo);
                
            case 'CHR (with 0% Overshoot)'
                [K, details] = designCHR(G, structure, options.epsilon, plantInfo);
                
            case 'Cohen-Coon'
                [K, details] = designCohenCoon(G, structure, options.epsilon, plantInfo);
                
            case 'Loop-Shaping'
                [K, details] = designLoopShaping(G, structure, options.phaseMargin, options.bandwidth, options.epsilon, plantInfo);
                
            case 'IMC (Internal Model Control)'
                [K, details] = designIMC(G, structure, options.settlingTime, options.epsilon, plantInfo);
                
            case 'MIGO (M-constrained Integral Gain Optimization)'
                [K, details] = designMIGO(G, structure, options.robustness, options.epsilon, plantInfo);
                
            case 'H-infinity'
                [K, details] = designHInfinity(G, structure, options.robustness, options.epsilon, plantInfo);
                
            case 'LQG (Linear-Quadratic-Gaussian)'
                [K, details] = designLQG(G, structure, options.bandwidth, options.robustness, options.epsilon, plantInfo);
                
            case 'Enhanced State Feedback'
                [K, details] = designEnhancedStateFeedback(G, structure, options, plantInfo);
                
            case 'Pre-stabilization'
                [K, details] = designPreStabilization(G, structure, options, plantInfo);
                
            case 'Youla-Kucera Parameterization'
                options.objective = options.goal;
                [K, details] = designYoulaKucera(G, structure, options, plantInfo);
                
            case 'Robust µ-synthesis'
                [K, details] = designMuSynthesis(G, structure, options, plantInfo);
                
            case 'Pole Placement'
                [K, details] = designPolePlacement(G, structure, options, plantInfo);
                
            otherwise
                error('Unknown design method: %s', method);
        end
    catch ME
        % Enhanced error handling with fallback controller design
        warning('Error in %s design: %s\nUsing robust fallback design.', method, ME.message);
        disp(['Design method error: ' ME.message]);
        [K, details] = designFallbackController(G, structure, options.epsilon, plantInfo);
        details = sprintf('Original design failed: %s\n\n%s', ME.message, details);
    end
    
    % Verify controller stability
    try
        [num, den] = tfdata(K, 'v');
        K_poles = roots(den);
        if any(real(K_poles) > 0)
            warning('Controller has unstable poles. Applying stabilization.');
            details = [details, '\nWARNING: Controller has unstable poles. Applying stabilization.\n'];
            
            % Stabilize the controller by reflecting unstable poles to LHP
            for i = 1:length(K_poles)
                if real(K_poles(i)) > 0
                    K_poles(i) = -abs(real(K_poles(i))) + imag(K_poles(i))*1i;
                end
            end
            
            % Create new controller with stabilized poles
            den_stable = poly(K_poles);
            K = tf(num, den_stable);
            details = [details, 'Controller poles have been stabilized.\n'];
        end
    catch
        % Continue if verification fails
    end
    
    % Verify closed-loop stability and adjust if needed
    try
        closed_loop = feedback(G*K, 1);
        cl_poles = pole(closed_loop);
        is_stable = all(real(cl_poles) < 0);
        
        if ~is_stable
            warning('Closed-loop system is unstable. Attempting gain adjustment.');
            details = [details, '\nWARNING: Closed-loop system is unstable. Attempting gain adjustment.\n'];
            
            % Try reducing gain until stable
            [num, den] = tfdata(K, 'v');
            stabilized = false;
            
            for scale = [0.5, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001]
                K_test = tf(num * scale, den);
                cl_test = feedback(G * K_test, 1);
                
                if all(real(pole(cl_test)) < 0)
                    K = K_test;
                    stabilized = true;
                    details = [details, sprintf('Closed-loop system stabilized with gain reduction factor: %.4f\n', scale)];
                    break;
                end
            end
            
            if ~stabilized
                details = [details, 'WARNING: Could not stabilize closed-loop system by gain reduction.\n'];
                details = [details, 'Trying fallback controller design...\n'];
                
                % If gain reduction fails, try more conservative controller
                [K_fallback, fallback_details] = designFallbackController(G, structure, options.epsilon, plantInfo);
                
                % Test fallback controller
                cl_fallback = feedback(G * K_fallback, 1);
                if all(real(pole(cl_fallback)) < 0)
                    K = K_fallback;
                    details = [details, 'FALLBACK CONTROLLER SUCCEEDED.\n\n'];
                    details = [details, fallback_details];
                else
                    details = [details, 'WARNING: Even fallback controller could not stabilize the system.\n'];
                    details = [details, 'Manual controller tuning is required for this challenging system.\n'];
                end
            end
        end
    catch ME
        details = [details, sprintf('\nError in closed-loop stability check: %s\n', ME.message)];
        details = [details, 'Stability could not be verified.\n'];
    end
    
    % Evaluate controller quality if requested
    try
        score = evaluateController(K, G, options.goal, options.phaseMargin, options.overshoot, options.settlingTime, options.bandwidth);
        details = [details, sprintf('\n\nController Score: %.2f/100', score)];
    catch ME
        disp(['Warning: Could not evaluate controller quality: ', ME.message]);
    end
end