function [K, details, score] = design_controller_auto(G, method, structure, options)
% DESIGN_CONTROLLER_AUTO Automatic controller design with various methods
% Enhanced version with intelligent method selection, state space integration,
% and robust error handling for highly unstable systems
% 
% Inputs:
%   G         - Plant model (transfer function or state-space)
%   method    - Design method (string) or 'auto' for automatic selection
%   structure - Controller structure: 'P', 'PI', 'PD', 'PID'
%   options   - Structure with optional parameters:
%     .epsilon     - Filter parameter for D-term (default: 0.1)
%     .phaseMargin - Desired phase margin in degrees (default: 45)
%     .bandwidth   - Desired bandwidth in rad/s (default: 1)
%     .settlingTime - Desired settling time in s (default: 5)
%     .robustness  - Robustness level: 'Low', 'Medium', 'High' (default: 'Medium')
%     .damping     - Desired damping ratio (default: 0.8)
%     .overshoot   - Desired overshoot in % (default: 10)
%     .goal        - Optimization goal: 'Tracking', 'Disturbance Rejection', 'Robustness' (default: 'Tracking')
%     .uncertainty - Uncertainty percentage for robust methods (default: 20)
%     .userSetBandwidth - Flag indicating if bandwidth was explicitly set by user (default: false)
%     .stateSpace  - Optional state-space model (if available from app)
%
% Outputs:
%   K        - Designed controller as transfer function
%   details  - Text description with design details
%   score    - Controller performance score (0-100)

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
        options.userSetBandwidth = false;
    elseif ~isfield(options, 'userSetBandwidth')
        options.userSetBandwidth = true;
    end
    
    if ~isfield(options, 'settlingTime')
        options.settlingTime = 5;
    end
    
    if ~isfield(options, 'robustness')
        options.robustness = 'Medium';
    end
    
    if ~isfield(options, 'damping')
        options.damping = 0.8;
    end
    
    if ~isfield(options, 'overshoot')
        options.overshoot = 10;
    end
    
    if ~isfield(options, 'goal')
        options.goal = 'Tracking';
    end
    
    if ~isfield(options, 'uncertainty')
        options.uncertainty = 20;
    end
    
    % Check if input is state-space or transfer function
    isStateSpace = isa(G, 'ss');
    
    % Store raw state-space model if available for methods that can use it directly
    if isStateSpace
        options.stateSpace = G;  % Store original state-space model
        [A, B, C, D] = ssdata(G);
        options.stateMatrices = {A, B, C, D}; % Store matrices separately
    end
    
    % Pre-analyze the plant to determine characteristics
    try
        plantInfo = analyzePlant(G);
        plantInfoString = getPlantInfoString(plantInfo);
        disp(['Plant Analysis: ', plantInfoString]);
    catch ME
        warning('Error in plant analysis: %s. Creating basic plant info.', ME.message);
        % Create basic plantInfo with default values
        plantInfo = createBasicPlantInfo(G);
    end
    
    % Check for highly unstable high-order systems
    if plantInfo.isUnstable && plantInfo.isHighOrder
        disp('Highly unstable high-order system detected. Adjusting design approach.');
        
        % Adjust parameters for such systems
        if ~options.userSetBandwidth
            % Calculate a safer bandwidth based on unstable poles
            p = plantInfo.poles;
            unstable_poles = p(real(p) > 0);
            max_real_part = max(real(unstable_poles));
            
            % Safer bandwidth value
            safe_bandwidth = max(0.2, min(options.bandwidth, max_real_part * 0.3));
            options.bandwidth = safe_bandwidth;
            disp(['Adjusted bandwidth to ', num2str(safe_bandwidth), ' rad/s for unstable system']);
        end
        
        % Increase robustness for unstable systems
        options.robustness = 'High';
        
        % Increase damping for more robust control
        options.damping = min(1.5, options.damping * 1.2);
        
        % Prioritize methods that work well with unstable systems
        if strcmpi(method, 'auto')
            disp('Automatically prioritizing methods suitable for unstable high-order systems');
            unstablePriorityMethods = {'Pre-stabilization', 'Enhanced State Feedback', 'Youla-Kucera Parameterization', 'Robust µ-synthesis'};
            method = selectBestMethodForUnstable(G, plantInfo, structure, unstablePriorityMethods);
        end
    end
    
    % Automatic method selection if method is 'auto'
    if strcmpi(method, 'auto')
        method = selectBestMethod(G, plantInfo, structure, options);
        disp(['Automatically selected method: ', method]);
    end
    
    % Prepare state space model if needed for specific methods
    stateSpaceNeededMethods = {'Enhanced State Feedback', 'Pre-stabilization', 'Pole Placement', 'LQG (Linear-Quadratic-Gaussian)'};
    
    if any(strcmpi(method, stateSpaceNeededMethods)) && ~isStateSpace && ~isfield(options, 'stateSpace')
        % Try to create state space model for methods that benefit from it
        try
            disp('Creating state-space model for advanced controller design...');
            options.stateSpace = ss(G);
            [A, B, C, D] = ssdata(options.stateSpace);
            options.stateMatrices = {A, B, C, D};
        catch ME
            disp(['Could not create standard state-space model: ', ME.message]);
            disp('Attempting enhanced state-space conversion...');
            
            try
                % Try enhanced numerical approach
                options.stateSpace = getEnhancedStateSpace(G, plantInfo);
                [A, B, C, D] = ssdata(options.stateSpace);
                options.stateMatrices = {A, B, C, D};
                disp('Successfully created enhanced state-space model');
            catch ME2
                disp(['Enhanced conversion also failed: ', ME2.message]);
                disp('Will continue without state-space model - may affect results');
            end
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
        
        % Special handling for unstable high-order systems
        if plantInfo.isUnstable && plantInfo.isHighOrder
            disp('Highly unstable system detected. Using specialized emergency approach.');
            [K, details] = designEmergencyControllerForUnstable(G, structure, options, plantInfo);
        else
            [K, details] = designFallbackController(G, structure, options.epsilon, plantInfo);
        end
        
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
    catch ME
        warning('Error while checking controller stability: %s', ME.message);
        details = [details, sprintf('\nWarning: Could not verify controller stability: %s\n', ME.message)];
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
                details = [details, 'Trying more advanced stabilization techniques...\n'];
                
                % Try more sophisticated approaches for unstable systems
                if plantInfo.isUnstable
                    [K_emergency, emergency_details] = designEmergencyControllerForUnstable(G, structure, options, plantInfo);
                    
                    % Test emergency controller
                    cl_emergency = feedback(G * K_emergency, 1);
                    if all(real(pole(cl_emergency)) < 0)
                        K = K_emergency;
                        details = [details, 'EMERGENCY CONTROLLER SUCCEEDED.\n\n'];
                        details = [details, emergency_details];
                    else
                        details = [details, 'WARNING: Even emergency controller could not stabilize the system.\n'];
                        details = [details, 'Manual controller tuning is required for this challenging system.\n'];
                    end
                else
                    % If system is originally stable, use very conservative controller
                    [K_fallback, fallback_details] = designConservativeController(G, structure, options.epsilon, plantInfo);
                    
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
        score = NaN;
    end
end

function method = selectBestMethod(G, plantInfo, structure, options)
    % SELECTBESTMETHOD Intelligently select the best control design method
    % based on plant characteristics and design requirements
    
    % Calculate a score for each method based on plant characteristics
    methods = {'Ziegler-Nichols (Oscillation)', 'Ziegler-Nichols (Step)', ...
              'Aström', 'CHR (with 0% Overshoot)', 'Cohen-Coon', ...
              'Loop-Shaping', 'IMC (Internal Model Control)', ...
              'MIGO (M-constrained Integral Gain Optimization)', ...
              'H-infinity', 'LQG (Linear-Quadratic-Gaussian)', ...
              'Enhanced State Feedback', 'Pre-stabilization', ...
              'Youla-Kucera Parameterization', 'Robust µ-synthesis', ...
              'Pole Placement'};
          
    % Initialize scores
    scores = zeros(1, length(methods));
    
    % Evaluate each method's suitability
    for i = 1:length(methods)
        scores(i) = evaluateMethodSuitability(methods{i}, G, plantInfo, structure, options);
    end
    
    % For unstable high-order systems, strongly bias toward specific methods
    if plantInfo.isUnstable && plantInfo.isHighOrder
        % These methods work best for challenging unstable systems
        unstable_methods = {'Pre-stabilization', 'Enhanced State Feedback', 'Youla-Kucera Parameterization', 'Robust µ-synthesis'};
        
        for i = 1:length(methods)
            if any(strcmpi(methods{i}, unstable_methods))
                scores(i) = scores(i) * 1.5;  % Boost score by 50%
            end
        end
    end
    
    % Find the best method (highest score)
    [~, bestIdx] = max(scores);
    method = methods{bestIdx};
    
    % Show top 3 methods and their scores for transparency
    [sortedScores, sortIdx] = sort(scores, 'descend');
    disp('Top 3 recommended methods:');
    for i = 1:min(3, length(methods))
        disp(sprintf('%d. %s (Score: %.2f)', i, methods{sortIdx(i)}, sortedScores(i)));
    end
end

function method = selectBestMethodForUnstable(G, plantInfo, structure, priorityMethods)
    % Special method selection for highly unstable systems
    
    % Test each priority method to see if it can be used
    for i = 1:length(priorityMethods)
        method = priorityMethods{i};
        
        try
            % Do a quick check to see if the method can handle this system
            if isApplicableForUnstable(G, method, plantInfo)
                disp(['Selected ' method ' as most suitable for this unstable system']);
                return;
            end
        catch
            % If method fails in preliminary check, try next one
            continue;
        end
    end
    
    % Default to pre-stabilization if all checks fail
    method = 'Pre-stabilization';
    disp('Defaulting to Pre-stabilization method as fallback for unstable system');
end

function applicable = isApplicableForUnstable(G, methodName, plantInfo)
    % Quick check if a method is applicable for unstable systems
    
    % Check method-specific criteria
    switch methodName
        case 'Pre-stabilization'
            % Always try pre-stabilization first for unstable systems
            applicable = true;
            
        case 'Enhanced State Feedback'
            % Check if state-space conversion is possible
            try
                ss(G);
                applicable = true;
            catch
                % If conversion fails, try enhanced conversion
                try
                    getEnhancedStateSpace(G, plantInfo);
                    applicable = true;
                catch
                    applicable = false;
                end
            end
            
        case 'Youla-Kucera Parameterization'
            % Check for some basic criteria
            p = plantInfo.poles;
            unstable_poles = p(real(p) > 0);
            
            % Youla-Kucera often struggles with too many unstable poles
            applicable = (length(unstable_poles) <= 3);
            
        case 'Robust µ-synthesis'
            % µ-synthesis may be too complex for very high order systems
            applicable = (length(plantInfo.poles) <= 10);
            
        otherwise
            applicable = false;
    end
end

function score = evaluateMethodSuitability(method, G, plantInfo, structure, options)
    % EVALUATEMETHODSUITABILITY Evaluate the suitability of a design method
    % for the given plant and requirements
    
    % Base score
    score = 50;
    
    % Plant characteristics-based scoring
    switch method
        case {'Ziegler-Nichols (Oscillation)', 'Ziegler-Nichols (Step)', 'Aström', 'CHR (with 0% Overshoot)', 'Cohen-Coon'}
            % Classical methods are good for stable, low-order systems
            if plantInfo.isUnstable
                score = score - 40;  % Greatly penalize for unstable plants
            end
            if plantInfo.hasRHPZeros
                score = score - 20;  % Penalize for non-minimum phase
            end
            if plantInfo.isHighOrder
                score = score - 15;  % Penalize for high-order systems
            end
            if plantInfo.hasIntegrator
                score = score - 10;  % Penalize for integrators (except ZN step)
                if strcmpi(method, 'Ziegler-Nichols (Step)')
                    score = score - 20;  % Further penalize ZN Step for integrators
                end
            end
            
        case 'Loop-Shaping'
            % Good general method, slightly penalize for complex systems
            if plantInfo.isUnstable
                score = score - 10;
            end
            if plantInfo.hasRHPZeros && plantInfo.isUnstable
                score = score - 20;  % Difficult combination
            end
            
        case 'IMC (Internal Model Control)'
            % Good for stable systems, especially with delays
            if plantInfo.isUnstable
                score = score - 30;
            end
            if plantInfo.hasDelay
                score = score + 15;  % Good for delay systems
            end
            
        case 'MIGO (M-constrained Integral Gain Optimization)'
            % Good for many systems requiring robustness
            if plantInfo.isUnstable && plantInfo.hasRHPZeros
                score = score - 20;
            end
            if strcmpi(options.robustness, 'High')
                score = score + 10;
            end
            
        case 'H-infinity'
            % Very good for robust control, especially with uncertainties
            if strcmpi(options.robustness, 'High')
                score = score + 15;
            end
            if plantInfo.hasRHPZeros
                score = score + 10;  % Good for non-minimum phase
            end
            if plantInfo.isUnstable
                score = score + 5;   % Can handle unstable plants
            end
            
        case 'LQG (Linear-Quadratic-Gaussian)'
            % Good for systems with noise, needs state-space representation
            if isa(G, 'ss') || isfield(options, 'stateSpace')
                score = score + 10;  % Bonus for state-space models
            end
            if plantInfo.isUnstable
                score = score + 5;   % Can handle some unstable plants
            end
            
        case 'Enhanced State Feedback'
            % Excellent for state-space models and unstable systems
            if isa(G, 'ss') || isfield(options, 'stateSpace')
                score = score + 15;
            end
            if plantInfo.isUnstable
                score = score + 25;  % Very good for unstable systems
            end
            if plantInfo.isHighOrder
                score = score + 15;  % Good for high-order systems
            end
            
        case 'Pre-stabilization'
            % Specifically designed for unstable systems
            if plantInfo.isUnstable
                score = score + 35;  % Major bonus for unstable systems
            else
                score = score - 25;  % Major penalty for stable systems
            end
            if plantInfo.hasRHPZeros && plantInfo.isUnstable
                score = score + 15;  % Extra bonus for this challenging combination
            end
            
        case 'Youla-Kucera Parameterization'
            % Great for complex and challenging systems
            if plantInfo.isUnstable || plantInfo.hasRHPZeros
                score = score + 25;
            end
            if plantInfo.isHighOrder
                score = score + 15;
            end
            
        case 'Robust µ-synthesis'
            % Excellent for dealing with uncertainties
            if strcmpi(options.robustness, 'High')
                score = score + 20;
            end
            if options.uncertainty > 15
                score = score + options.uncertainty / 2;  % Scales with uncertainty level
            end
            if plantInfo.isUnstable
                score = score + 15;  % Good for unstable systems
            end
            
        case 'Pole Placement'
            % Good for state-space models and when precise dynamics are needed
            if isa(G, 'ss') || isfield(options, 'stateSpace')
                score = score + 15;
            end
            if plantInfo.isUnstable
                score = score + 20;  % Very good for unstable systems
            end
            if plantInfo.isHighOrder
                score = score + 5;   % Can handle high-order systems
            end
    end
    
    % Controller structure-based adjustments
    if strcmpi(structure, 'P')
        % Simpler methods may be better for P controllers
        if ismember(method, {'Ziegler-Nichols (Oscillation)', 'Loop-Shaping', 'MIGO (M-constrained Integral Gain Optimization)'})
            score = score + 5;
        end
      elseif strcmpi(structure, 'PID')
        % More advanced methods may be better for PID controllers
        if ismember(method, {'Youla-Kucera Parameterization', 'H-infinity', 'Robust µ-synthesis', 'Enhanced State Feedback'})
            score = score + 5;
        end
    end
    
    % Performance requirements-based adjustments
    if strcmpi(options.goal, 'Tracking')
        if ismember(method, {'IMC (Internal Model Control)', 'Pole Placement', 'Enhanced State Feedback'})
            score = score + 5;
        end
      elseif strcmpi(options.goal, 'Disturbance Rejection')
        if ismember(method, {'MIGO (M-constrained Integral Gain Optimization)', 'H-infinity', 'Youla-Kucera Parameterization'})
            score = score + 5;
        end
      elseif strcmpi(options.goal, 'Robustness')
        if ismember(method, {'H-infinity', 'Robust µ-synthesis', 'Youla-Kucera Parameterization'})
            score = score + 10;
        end
    end
    
    % Bandwidth considerations
    if options.bandwidth > 5 && ismember(method, {'Pre-stabilization', 'Enhanced State Feedback', 'Robust µ-synthesis'})
        score = score + 5;  % These methods can handle high-bandwidth requirements
    end
    
    % Balance complexity vs. effectiveness for unstable systems
    if plantInfo.isUnstable && length(plantInfo.poles) > 4
        if ismember(method, {'Pole Placement', 'Enhanced State Feedback', 'Pre-stabilization'})
            score = score + 15;  % These methods work well for high-order unstable systems
        end
    end
    
    % Final score adjustments
    score = min(max(score, 0), 100);  % Ensure score is between 0 and 100
end

function plantInfo = createBasicPlantInfo(G)
    % CREATEBASICPLANTINFO Create a basic plantInfo structure when
    % analyzePlant fails
    
    plantInfo = struct();
    
    % Try to get basic plant properties
    try
        plantInfo.poles = pole(G);
    catch
        plantInfo.poles = [0];  % Default if poles can't be computed
    end
    
    try
        plantInfo.zeros = zero(G);
    catch
        plantInfo.zeros = [];  % Default if zeros can't be computed
    end
    
    % Determine basic stability and properties
    plantInfo.isUnstable = any(real(plantInfo.poles) > 0);
    plantInfo.hasIntegrator = any(abs(plantInfo.poles) < 1e-6);
    plantInfo.hasRHPZeros = any(real(plantInfo.zeros) > 0);
    plantInfo.hasDelay = false;  % Default assumption
    plantInfo.isHighOrder = length(plantInfo.poles) > 2;
    
    % Try to get DC gain
    try
        plantInfo.dcGain = dcgain(G);
    catch
        plantInfo.dcGain = NaN;
    end
    
    % Add default values for other fields
    plantInfo.stepInfo = [];
    plantInfo.stepResponse = [];
    plantInfo.FOPDT = struct('L', NaN, 'T', NaN, 'K', NaN);
    plantInfo.stabilizingGain = 0;
    if plantInfo.isUnstable
        % Estimate stabilizing gain for unstable plants
        try
            plantInfo.stabilizingGain = estimateStabilizingGain(G);
        catch
            plantInfo.stabilizingGain = 1;  % Default if estimation fails
        end
    end
    
    % Add empty controls for other fields from analyzePlant
    plantInfo.controllability = struct('rank', NaN, 'size', NaN, 'condition', NaN);
    plantInfo.observability = struct('rank', NaN, 'size', NaN, 'condition', NaN);
    plantInfo.modes = struct('eigenvalues', [], 'controllability', [], 'observability', []);
    plantInfo.poorlyControllableModes = [];
    plantInfo.poorlyObservableModes = [];
end

function infoStr = getPlantInfoString(plantInfo)
    % Create a formatted string with plant information
    
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

function [K, details] = designEmergencyControllerForUnstable(G, structure, options, plantInfo)
    % Design emergency controller for highly unstable systems when all else fails
    
    details = 'EMERGENCY CONTROLLER DESIGN FOR HIGHLY UNSTABLE SYSTEM\n';
    details = [details, '---------------------------------------------------\n'];
    details = [details, sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo))];
    
    % First try direct pre-stabilization as pure pole-zero cancellation
    try
        % Extract poles for direct cancellation
        p = plantInfo.poles;
        unstable_poles = p(real(p) > 0);
        
        K_num = 1;
        K_den = 1;
        
        % Direct pole-zero cancellation
        for i = 1:length(unstable_poles)
            pole_i = unstable_poles(i);
            
            if imag(pole_i) ~= 0
                % Skip conjugate pairs, we'll add both together
                if i < length(unstable_poles) && abs(pole_i - conj(unstable_poles(i+1))) < 1e-6
                    continue;
                end
                
                if imag(pole_i) > 0
                    real_part = real(pole_i);
                    imag_part = imag(pole_i);
                    
                    % Create zeros at unstable poles
                    quad_term = [1, -2*real_part, real_part^2 + imag_part^2];
                    
                    % Create stable poles at reflected locations with extra damping
                    stable_quad = [1, 5*real_part, 6*(real_part^2 + imag_part^2)];
                    
                    K_num = conv(K_num, quad_term);
                    K_den = conv(K_den, stable_quad);
                end
              else
                % For real poles
                K_num = conv(K_num, [1, -pole_i]);
                K_den = conv(K_den, [1, 5*pole_i]);  % More damping for stability
            end
        end
        
        % Use extremely conservative gain for highly unstable systems
        gain_factor = 0.001;
        
        K_stab = tf(gain_factor * K_num, K_den);
        
        % Verify stabilization
        closed_loop = feedback(G * K_stab, 1);
        cl_poles = pole(closed_loop);
        
        if all(real(cl_poles) < 0)
            details = [details, 'Direct pole cancellation successful!\n'];
            
            % For this emergency case, no need for advanced controller structure
            % Just use the stabilizing controller
            K = K_stab;
            details = [details, 'Using stabilizing controller directly due to high instability.\n'];
            return;
          else
            details = [details, 'Direct cancellation failed. Trying gain adjustments...\n'];
            
            % Try different gain scaling factors
            for scale = [0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001]
                K_test = tf(gain_factor * scale * K_num, K_den);
                closed_loop = feedback(G * K_test, 1);
                
                if all(real(pole(closed_loop)) < 0)
                    K = K_test;
                    details = [details, sprintf('System stabilized with gain factor %.6f\n', gain_factor * scale)];
                    return;
                end
            end
        end
      catch ME
        details = [details, sprintf('Pole cancellation approach failed: %s\n', ME.message)];
    end
    
    % If direct cancellation fails, try with state space approach
    try
        details = [details, 'Attempting state-space stabilization...\n'];
        
        % Try to get or create state-space model
        if isfield(options, 'stateSpace') && ~isempty(options.stateSpace)
            sys_ss = options.stateSpace;
          else
            sys_ss = getEnhancedStateSpace(G, plantInfo);
        end
        
        % Extract state-space matrices
        A = sys_ss.A;
        B = sys_ss.B;
        C = sys_ss.C;
        n = size(A, 1);
        
        % Calculate very conservative poles far in left half-plane
        p = eig(A);
        unstable_idx = find(real(p) > 0);
        
        if ~isempty(unstable_idx)
            % Create desired_poles by reflecting unstable poles
            desired_poles = p;
            
            for i = 1:length(unstable_idx)
                idx = unstable_idx(i);
                if imag(p(idx)) ~= 0
                    % For complex poles, maintain frequency but add heavy damping
                    freq = abs(p(idx));
                    desired_poles(idx) = -freq * 0.8 + 1j * freq * 0.2;  % 80% damping
                  else
                    % For real poles, reflect with margin
                    desired_poles(idx) = -abs(real(p(idx))) * 5;  % 5x margin
                end
            end
            
            % Make sure we have conjugate pairs
            for i = 1:length(desired_poles)
                if imag(desired_poles(i)) ~= 0
                    % Find if this pole has a conjugate pair
                    has_conjugate = false;
                    for j = 1:length(desired_poles)
                        if i ~= j && abs(desired_poles(i) - conj(desired_poles(j))) < 1e-10
                            has_conjugate = true;
                            break;
                        end
                    end
                    
                    % If no conjugate found, replace a real pole with the conjugate
                    if ~has_conjugate
                        for j = 1:length(desired_poles)
                            if imag(desired_poles(j)) == 0
                                desired_poles(j) = conj(desired_poles(i));
                                break;
                            end
                        end
                    end
                end
            end
            
            % Try different pole placement methods
            try
                K_state = place(A, B, desired_poles);
              catch
                try
                    K_state = acker(A, B, desired_poles);
                  catch
                    % Last resort: simple approximation
                    K_state = zeros(1, n);
                    for i = 1:length(unstable_idx)
                        idx = unstable_idx(i);
                        pole_i = p(idx);
                        shift_vector = zeros(n, 1);
                        shift_vector(idx) = 1;
                        shift_amount = real(pole_i) * 10;  % 10x margin
                        K_state = K_state + shift_amount * shift_vector' * pinv(B);
                    end
                end
            end
            
            % Create state feedback controller
            Ac = A - B*K_state;
            Bc = B;
            Cc = -K_state;
            Dc = 0;
            
            K_ss = ss(Ac, Bc, Cc, Dc);
            K = tf(K_ss);
            
            % Verify stabilization
            closed_loop = feedback(G*K, 1);
            cl_poles = pole(closed_loop);
            
            if all(real(cl_poles) < 0)
                details = [details, 'State-space stabilization successful!\n'];
                return;
              else 
                details = [details, 'State-space approach failed. Trying gain adjustments...\n'];
                
                % Try scaling the gain
                [num, den] = tfdata(K, 'v');
                
                for scale = [0.5, 0.1, 0.05, 0.01, 0.005, 0.001]
                    K_test = tf(num * scale, den);
                    closed_loop = feedback(G * K_test, 1);
                    
                    if all(real(pole(closed_loop)) < 0)
                        K = K_test;
                        details = [details, sprintf('System stabilized with gain scale factor %.4f\n', scale)];
                        return;
                    end
                end
            end
        end
      catch ME
        details = [details, sprintf('State-space approach failed: %s\n', ME.message)];
    end
    
    % Last resort: create a very basic controller based on structure
    details = [details, 'All approaches failed. Creating ultra-conservative controller.\n'];
    
    switch structure
        case 'P'
            K = tf(0.001, 1);
        case 'PI'
            K = tf([0.001, 0.0001], [1, 0]);
        case 'PD'
            K = tf([0.01, 0.001], [0.001, 1]);
        case 'PID'
            K = tf([0.01, 0.001, 0.0001], [0.001, 1, 0]);
        otherwise
            K = tf(0.001, 1);
    end
    
    details = [details, 'WARNING: Ultra-conservative controller created. Manual tuning strongly recommended.\n'];
    return;
end

function [K, details] = designConservativeController(G, structure, epsilon, plantInfo)
    % Design very conservative controller when other methods fail
    
    details = 'CONSERVATIVE CONTROLLER DESIGN\n';
    details = [details, '----------------------------\n'];
    details = [details, sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo))];
    
    % Get a very safe gain value
    try
        dc_gain = abs(dcgain(G));
        if isnan(dc_gain) || isinf(dc_gain) || dc_gain == 0
            dc_gain = 1;
        end
        
        safe_gain = 0.1 / dc_gain;
      catch
        safe_gain = 0.01;
    end
    
    % Limit to a reasonable range
    safe_gain = min(max(safe_gain, 0.001), 1);
    
    % Create simple controllers with very conservative parameters
    switch structure
        case 'P'
            K = tf(safe_gain, 1);
            details = [details, sprintf('Created conservative P controller: Kp = %.6f\n', safe_gain)];
            
        case 'PI'
            Kp = safe_gain;
            Ki = safe_gain * 0.05;  % Very low integral action
            K = tf([Kp, Ki], [1, 0]);
            details = [details, sprintf('Created conservative PI controller: Kp = %.6f, Ki = %.6f\n', Kp, Ki)];
            
        case 'PD'
            Kp = safe_gain;
            Kd = safe_gain * 0.5;
            K = tf([Kd, Kp], [epsilon*10, 1]);  % Extra filtering
            details = [details, sprintf('Created conservative PD controller: Kp = %.6f, Kd = %.6f\n', Kp, Kd)];
            
        case 'PID'
            Kp = safe_gain;
            Ki = safe_gain * 0.01;  % Very low integral action
            Kd = safe_gain * 0.5;
            K = tf([Kd, Kp, Ki], [epsilon*10, 1, 0]);  % Extra filtering
            details = [details, sprintf('Created conservative PID controller: Kp = %.6f, Ki = %.6f, Kd = %.6f\n', Kp, Ki, Kd)];
            
        otherwise
            K = tf(safe_gain, 1);
            details = [details, sprintf('Created conservative controller with gain = %.6f\n', safe_gain)];
    end
    
    % Verify stability
    try
        closed_loop = feedback(G*K, 1);
        cl_poles = pole(closed_loop);
        
        is_stable = all(real(cl_poles) < 0);
        
        if is_stable
            details = [details, 'Conservative controller successfully stabilizes the system.\n'];
          else 
            details = [details, 'WARNING: Even conservative controller does not stabilize the system.\n'];
            details = [details, 'This system is extremely challenging to control.\n'];
            
            % Try reducing gain even further
            for scale = [0.1, 0.01, 0.001, 0.0001]
                [num, den] = tfdata(K, 'v');
                K_test = tf(num * scale, den);
                closed_loop = feedback(G * K_test, 1);
                
                if all(real(pole(closed_loop)) < 0)
                    K = K_test;
                    details = [details, sprintf('System stabilized with additional gain reduction (factor: %.5f)\n', scale)];
                    break;
                end
            end
        end
      catch ME
        details = [details, sprintf('Could not verify stability: %s\n', ME.message)];
    end
    
    return;
end