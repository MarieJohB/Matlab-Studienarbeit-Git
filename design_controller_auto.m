function [K, details, score] = design_controller_auto(G, method, structure, options)
% DESIGN_CONTROLLER_AUTO Automatic controller design with various methods
% Enhanced version with intelligent method selection and robust error handling
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
    
    % Automatic method selection if method is 'auto'
    if strcmpi(method, 'auto')
        method = selectBestMethod(G, plantInfo, structure, options);
        disp(['Automatically selected method: ', method]);
    end
    
    % For highly unstable systems, adjust parameters
    if plantInfo.isUnstable
        % Get a measure of instability
        unstable_poles = plantInfo.poles(real(plantInfo.poles) > 0);
        max_real_part = max(real(unstable_poles));
        
        % For highly unstable systems, use more conservative settings
        if max_real_part > 5 || length(unstable_poles) > 1
            disp('System is highly unstable - adjusting control parameters for better stability');
            
            % Reduce bandwidth for very unstable systems 
            if ~options.userSetBandwidth
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
            
        case 'LQG (Linear-Quadratic-Gaussian)'
            % Good for systems with noise, needs state-space representation
            if isa(G, 'ss')
                score = score + 10;  % Bonus for state-space models
            end
            
        case 'Enhanced State Feedback'
            % Excellent for state-space models and unstable systems
            if isa(G, 'ss')
                score = score + 15;
            end
            if plantInfo.isUnstable
                score = score + 15;  % Very good for unstable systems
            end
            if plantInfo.isHighOrder
                score = score + 10;  % Good for high-order systems
            end
            
        case 'Pre-stabilization'
            % Specifically designed for unstable systems
            if plantInfo.isUnstable
                score = score + 25;  % Major bonus for unstable systems
            else
                score = score - 25;  % Major penalty for stable systems
            end
            if plantInfo.hasRHPZeros && plantInfo.isUnstable
                score = score + 15;  % Extra bonus for this challenging combination
            end
            
        case 'Youla-Kucera Parameterization'
            % Great for complex and challenging systems
            if plantInfo.isUnstable || plantInfo.hasRHPZeros
                score = score + 20;
            end
            if plantInfo.isHighOrder
                score = score + 10;
            end
            
        case 'Robust µ-synthesis'
            % Excellent for dealing with uncertainties
            if strcmpi(options.robustness, 'High')
                score = score + 20;
            end
            if options.uncertainty > 15
                score = score + options.uncertainty / 2;  % Scales with uncertainty level
            end
            
        case 'Pole Placement'
            % Good for state-space models and when precise dynamics are needed
            if isa(G, 'ss')
                score = score + 15;
            end
            if plantInfo.isUnstable
                score = score + 10;
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