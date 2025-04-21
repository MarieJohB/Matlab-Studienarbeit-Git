function [K, details, score] = design_controller_auto(G, method, structure, options)
% DESIGN_CONTROLLER_AUTO Automatic controller design with enhanced robustness
% Improved version with intelligent method selection, enhanced error handling, 
% and multi-level fallback mechanisms for highly unstable systems
% 
% Inputs:
%   G         - Plant model (transfer function or state-space)
%   method    - Design method or 'auto' for automatic selection
%   structure - Controller structure: 'P', 'PI', 'PD', 'PID'
%   options   - Structure with optional parameters
%
% Outputs:
%   K        - Designed controller as transfer function
%   details  - Text description with design details
%   score    - Controller performance score (0-100)

    % Default values for missing options
    if nargin < 4
        options = struct();
    end
    
    % Set default options with improved initialization
    options = setDefaultOptions(options);
    
    % Check if input is state-space or transfer function
    isStateSpace = isa(G, 'ss');
    
    % Store raw state-space model if available for methods that use it directly
    if isStateSpace
        options.stateSpace = G;  % Store original state-space model
        [A, B, C, D] = ssdata(G);
        options.stateMatrices = {A, B, C, D}; % Store matrices separately
        
        % Convert to transfer function for methods that require it
        try
            G_tf = tf(G);
        catch ME
            warning('State-space to transfer function conversion failed: %s', ME.message);
            G_tf = G;  % Keep as state-space and let methods handle it
        end
    else
        G_tf = G;  % Already a transfer function
    end
    
    % Pre-analyze the plant to determine characteristics
    try
        plantInfo = analyzePlant(G_tf);
        plantInfoString = getPlantInfoString(plantInfo);
        disp(['Plant Analysis: ', plantInfoString]);
    catch ME
        warning('Error in plant analysis: %s. Creating basic plant info.', ME.message);
        % Create basic plantInfo with default values
        plantInfo = createBasicPlantInfo(G_tf);
    end
    
    % Enhanced classification of difficult plants 
    [plantDifficulty, specialConditions] = classifyPlantDifficulty(plantInfo);
    disp(['Plant difficulty classification: ', plantDifficulty]);
    
    % If special conditions exist, display them
    if ~isempty(fieldnames(specialConditions))
        disp('Special plant conditions detected:');
        fields = fieldnames(specialConditions);
        for i = 1:length(fields)
            if specialConditions.(fields{i})
                disp(['  - ', fields{i}]);
            end
        end
    end
    
    % Adjust options based on plant characteristics for safer design
    % FIXED: Pass structure as parameter to adjustOptionsForPlantType
    options = adjustOptionsForPlantType(options, plantInfo, plantDifficulty, specialConditions, structure);
    
    % Automatic method selection if method is 'auto'
    if strcmpi(method, 'auto')
        method = selectBestMethod(G_tf, plantInfo, structure, options, plantDifficulty, specialConditions);
        disp(['Automatically selected method: ', method]);
    end
    
    % Prepare state space model if needed for specific methods
    stateSpaceNeededMethods = {'Pole Placement', 'LQR', 'LQG'};
    
    if any(strcmpi(method, stateSpaceNeededMethods)) && ~isStateSpace && ~isfield(options, 'stateSpace')
        options = prepareStateSpaceModel(G_tf, options, plantInfo);
    end
    
    % Initial design attempt with selected method
    try
        [K, details] = designControllerWithMethod(method, G_tf, structure, options, plantInfo);
        designSuccess = true;
    catch ME
        % First level fallback on method failure
        warning('Primary design method (%s) failed: %s\nAttempting fallback method.', method, ME.message);
        disp(['Design method error: ' ME.message]);
        designSuccess = false;
        
        % Store original error for details
        originalError = ME.message;
    end
    
    % If initial design failed, try fallback methods
    if ~designSuccess
        [K, details, fallbackSuccess] = attemptFallbackDesign(G_tf, structure, options, plantInfo, plantDifficulty, originalError);
        
        % If all fallbacks failed, create a minimal stabilizing controller
        if ~fallbackSuccess
            [K, details] = createEmergencyController(G_tf, structure, options, plantInfo, plantDifficulty);
        end
    end
    
    % Post-process controller to ensure stability and proper structure
    [K, details] = postProcessController(K, G_tf, details, structure, options, plantInfo);
    
    % Evaluate controller quality
    try
        score = evaluateController(K, G_tf, options.goal, options.phaseMargin, options.overshoot, options.settlingTime, options.bandwidth, plantInfo);
        details = [details, sprintf('\n\nController Score: %.2f/100', score)];
    catch ME
        disp(['Warning: Could not evaluate controller quality: ', ME.message]);
        score = NaN;
    end
end

%% Helper Functions

function options = setDefaultOptions(options)
    % Set default options with safe values
    
    % Essential parameters
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
    
    % Robustness and performance parameters
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
    
    % Ensure valid robustness level
    validRobustnessLevels = {'Low', 'Medium', 'High', 'Very High'};
    if ~ismember(options.robustness, validRobustnessLevels)
        options.robustness = 'Medium';
    end
end

function [difficulty, specialConditions] = classifyPlantDifficulty(plantInfo)
    % Classify plant difficulty based on detailed analysis
    
    % Initialize special conditions structure
    specialConditions = struct();
    specialConditions.multipleUnstablePoles = false;
    specialConditions.highlyUnstable = false;
    specialConditions.closeRHPZeros = false;
    specialConditions.complexRHPZeros = false;
    specialConditions.rhpZeroNearUnstablePole = false;
    specialConditions.highlyOscillatory = false;
    specialConditions.integrationWithRHP = false;
    specialConditions.numericallyIll = false;
    specialConditions.nonCollocated = false;
    
    % Get poles and zeros
    p = plantInfo.poles;
    z = plantInfo.zeros;
    
    % Check for unstable poles
    unstable_poles = p(real(p) > 0);
    specialConditions.multipleUnstablePoles = length(unstable_poles) > 1;
    
    % Check for highly unstable poles
    if ~isempty(unstable_poles)
        max_real_pole = max(real(unstable_poles));
        specialConditions.highlyUnstable = max_real_pole > 5;
    end
    
    % Check for RHP zeros
    rhp_zeros = z(real(z) > 0);
    
    if ~isempty(rhp_zeros)
        % Check for RHP zeros close to imaginary axis
        min_real_zero = min(real(rhp_zeros));
        specialConditions.closeRHPZeros = min_real_zero < 0.5;
        
        % Check for complex RHP zeros
        specialConditions.complexRHPZeros = any(imag(rhp_zeros) ~= 0);
        
        % Check for RHP zeros near unstable poles (very challenging condition)
        if ~isempty(unstable_poles)
            for i = 1:length(rhp_zeros)
                for j = 1:length(unstable_poles)
                    if abs(rhp_zeros(i) - unstable_poles(j)) < 2
                        specialConditions.rhpZeroNearUnstablePole = true;
                        break;
                    end
                end
                if specialConditions.rhpZeroNearUnstablePole
                    break;
                end
            end
        end
    end
    
    % Check for highly oscillatory dynamics
    oscillatory_poles = p(abs(imag(p)) > abs(real(p)));
    specialConditions.highlyOscillatory = ~isempty(oscillatory_poles);
    
    % Check for integration with RHP dynamics
    specialConditions.integrationWithRHP = plantInfo.hasIntegrator && (plantInfo.isUnstable || plantInfo.hasRHPZeros);
    
    % Check for numerical conditioning
    try
        [num, den] = tfdata(tf(plantInfo.poles, plantInfo.zeros, 1), 'v');
        max_coeff = max(max(abs(num)), max(abs(den)));
        min_num = min(abs(num(abs(num) > 0)));
        min_den = min(abs(den(abs(den) > 0)));
        min_coeff = min(min_num, min_den);
        
        condition_number = max_coeff / min_coeff;
        specialConditions.numericallyIll = condition_number > 1e8;
    catch
        % Default to false if calculation fails
    end
    
    % Determine overall difficulty
    if specialConditions.multipleUnstablePoles && (specialConditions.closeRHPZeros || specialConditions.rhpZeroNearUnstablePole)
        difficulty = 'Extremely Difficult';
    elseif specialConditions.multipleUnstablePoles || specialConditions.highlyUnstable
        difficulty = 'Very Difficult';
    elseif specialConditions.rhpZeroNearUnstablePole || (plantInfo.isUnstable && plantInfo.hasRHPZeros)
        difficulty = 'Very Difficult';
    elseif plantInfo.isUnstable || (specialConditions.closeRHPZeros && specialConditions.complexRHPZeros)
        difficulty = 'Difficult';
    elseif plantInfo.hasRHPZeros || specialConditions.highlyOscillatory || specialConditions.integrationWithRHP
        difficulty = 'Moderately Difficult';
    elseif plantInfo.isHighOrder || plantInfo.hasIntegrator || specialConditions.numericallyIll
        difficulty = 'Somewhat Challenging';
    else
        difficulty = 'Standard';
    end
end

function options = adjustOptionsForPlantType(options, plantInfo, plantDifficulty, specialConditions, structure)
    % ADJUSTOPTIONSFORPLANTTYPE Adjust design options based on plant characteristics
    % 
    % Inputs:
    %   options          - Design options structure
    %   plantInfo        - Plant information structure
    %   plantDifficulty  - Plant difficulty classification
    %   specialConditions - Special plant conditions structure
    %   structure        - Controller structure ('P', 'PI', 'PD', 'PID')
    %
    % Output:
    %   options - Adjusted design options
    
    % Set default options
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
    
    % Ensure valid robustness level
    validRobustnessLevels = {'Low', 'Medium', 'High', 'Very High'};
    if ~ismember(options.robustness, validRobustnessLevels)
        options.robustness = 'Medium';
    end
    
    % For very difficult plants, adjust bandwidth if user didn't explicitly set it
    % FIXED: Changed 'difficulty' to 'plantDifficulty'
    switch plantDifficulty
        case 'Extremely Difficult'
            % Very conservative bandwidth for extremely difficult plants
            if ~options.userSetBandwidth
                if specialConditions.multipleUnstablePoles
                    % For multiple unstable poles, calculate safer bandwidth
                    unstable_poles = plantInfo.poles(real(plantInfo.poles) > 0);
                    max_real_part = max(real(unstable_poles));
                    safe_bandwidth = max(0.1, min(options.bandwidth, max_real_part * 0.2));
                else
                    safe_bandwidth = min(options.bandwidth, 0.3);
                end
                options.bandwidth = safe_bandwidth;
                disp(['Adjusted bandwidth to ', num2str(safe_bandwidth), ' rad/s for extremely difficult plant']);
            end
            
            % Increase damping for stability
            options.damping = max(options.damping, 0.9);
            
            % Increase phase margin for stability
            options.phaseMargin = max(options.phaseMargin, 60);
            
            % Higher epsilon for derivative filtering
            if ismember(structure, {'PD', 'PID'})
                options.epsilon = max(options.epsilon, 0.2);
            end
            
        case 'Very Difficult'
            % Conservative bandwidth for very difficult plants
            if ~options.userSetBandwidth
                if plantInfo.isUnstable
                    unstable_poles = plantInfo.poles(real(plantInfo.poles) > 0);
                    max_real_part = max(real(unstable_poles));
                    safe_bandwidth = max(0.2, min(options.bandwidth, max_real_part * 0.3));
                else
                    safe_bandwidth = min(options.bandwidth, 0.5);
                end
                options.bandwidth = safe_bandwidth;
                disp(['Adjusted bandwidth to ', num2str(safe_bandwidth), ' rad/s for very difficult plant']);
            end
            
            % Increased damping
            options.damping = max(options.damping, 0.8);
            
            % Increased phase margin
            options.phaseMargin = max(options.phaseMargin, 55);
            
            % Higher epsilon for derivative filtering
            if ismember(structure, {'PD', 'PID'})
                options.epsilon = max(options.epsilon, 0.15);
            end
            
        case 'Difficult'
            % Somewhat conservative bandwidth for difficult plants
            if ~options.userSetBandwidth && plantInfo.hasRHPZeros
                rhp_zeros = plantInfo.zeros(real(plantInfo.zeros) > 0);
                min_real_part = min(real(rhp_zeros));
                safe_bandwidth = min(options.bandwidth, min_real_part * 0.4);
                options.bandwidth = safe_bandwidth;
                disp(['Adjusted bandwidth to ', num2str(safe_bandwidth), ' rad/s due to RHP zero limitations']);
            elseif ~options.userSetBandwidth && plantInfo.isUnstable
                unstable_poles = plantInfo.poles(real(plantInfo.poles) > 0);
                max_real_part = max(real(unstable_poles));
                safe_bandwidth = min(options.bandwidth, max_real_part * 0.5);
                options.bandwidth = safe_bandwidth;
                disp(['Adjusted bandwidth to ', num2str(safe_bandwidth), ' rad/s for unstable plant']);
            end
            
            % Moderately increased phase margin
            options.phaseMargin = max(options.phaseMargin, 50);
            
            % Moderate epsilon for derivative filtering
            if ismember(structure, {'PD', 'PID'})
                options.epsilon = max(options.epsilon, 0.1);
            end
    end
    
    % Adjust robustness based on difficulty
    % FIXED: Changed 'difficulty' to 'plantDifficulty'
    switch plantDifficulty
        case 'Extremely Difficult'
            options.robustness = 'Very High';
        case 'Very Difficult'
            options.robustness = 'High';
        case 'Difficult'
            if ~strcmpi(options.robustness, 'Very High')
                options.robustness = 'High';
            end
    end
end

function method = selectBestMethod(G, plantInfo, structure, options, difficulty, specialConditions)
    % Enhanced method selection with improved handling of difficult plants
    
    % Available design methods
    methods = {'Ziegler-Nichols (Oscillation)', 'Ziegler-Nichols (Step)', ...
              'Aström', 'Loop-Shaping', 'Pole Placement', 'Compensation Controller'};
    
    % Prefilter methods for extreme difficulty levels
    switch difficulty
        case 'Extremely Difficult'
            % For extremely difficult plants, prefer only specialized methods
            preferredMethods = {'Compensation Controller', 'Pole Placement'};
            disp('Limiting method selection to specialized techniques for extremely difficult plant');
            
            % Check if state-space representation is available for pole placement
            if isfield(options, 'stateSpace') || isa(G, 'ss')
                % Prefer state-space method if available
                method = 'Pole Placement';
            else
                % Default to compensation controller
                method = 'Compensation Controller';
            end
            return;
            
        case 'Very Difficult'
            % For very difficult plants, prioritize these methods
            preferredMethods = {'Compensation Controller', 'Pole Placement', 'Loop-Shaping'};
            
            % Special case: multiple unstable poles without state-space model
            if specialConditions.multipleUnstablePoles && ~(isfield(options, 'stateSpace') || isa(G, 'ss'))
                method = 'Compensation Controller';
                disp('Selected Compensation Controller for multiple unstable poles');
                return;
            end
    end
    
    % Evaluate method suitability for standard cases
    scores = zeros(1, length(methods));
    
    for i = 1:length(methods)
        scores(i) = evaluateMethodSuitability(methods{i}, G, plantInfo, structure, options, difficulty, specialConditions);
    end
    
    % Special handling for difficult plants
    if strcmp(difficulty, 'Very Difficult') || strcmp(difficulty, 'Difficult')
        % Boost scores for preferred methods for difficult plants
        for i = 1:length(methods)
            if ismember(methods{i}, preferredMethods)
                scores(i) = scores(i) * 1.5; % 50% boost
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

function score = evaluateMethodSuitability(method, G, plantInfo, structure, options, difficulty, specialConditions)
    % Evaluate how suitable a design method is for the given plant and requirements
    
    % Base score
    score = 50;
    
    % Enhanced plant characteristics-based scoring
    switch method
        case {'Ziegler-Nichols (Oscillation)', 'Ziegler-Nichols (Step)', 'Aström'}
            % Classical methods are good for simple systems
            if strcmp(difficulty, 'Extremely Difficult')
                score = score - 45;  % Practically eliminate for extremely difficult plants
            elseif strcmp(difficulty, 'Very Difficult')
                score = score - 40;  % Strongly penalize for very difficult plants
            elseif strcmp(difficulty, 'Difficult')
                score = score - 30;  % Heavily penalize for difficult plants
            end
            
            % Special conditions penalties
            if specialConditions.multipleUnstablePoles
                score = score - 40;  % Severe penalty for multiple unstable poles
            end
            if specialConditions.highlyUnstable
                score = score - 35;  % Strong penalty for highly unstable poles
            end
            if specialConditions.closeRHPZeros
                score = score - 25;  % Penalty for close RHP zeros
            end
            if specialConditions.rhpZeroNearUnstablePole
                score = score - 40;  % Severe penalty for RHP zeros near unstable poles
            end
            if specialConditions.integrationWithRHP
                score = score - 30;  % Strong penalty for integrators with RHP dynamics
            end
            
            % ZN Step specifically struggles with integrators
            if strcmpi(method, 'Ziegler-Nichols (Step)') && plantInfo.hasIntegrator
                score = score - 25;  % Extra penalty
            end
            
        case 'Loop-Shaping'
            % Loop shaping can handle a variety of systems but struggles with very challenging ones
            if strcmp(difficulty, 'Extremely Difficult')
                score = score - 35;  % Heavy penalty for extremely difficult plants
            elseif strcmp(difficulty, 'Very Difficult')
                score = score - 20;  % Moderate penalty for very difficult plants
            elseif strcmp(difficulty, 'Difficult')
                score = score - 10;  % Small penalty for difficult plants
            end
            
            % Special conditions penalties (and bonuses)
            if specialConditions.multipleUnstablePoles
                score = score - 30;  % Strong penalty for multiple unstable poles
            end
            if specialConditions.highlyUnstable
                score = score - 25;  % Penalty for highly unstable poles
            end
            if specialConditions.closeRHPZeros && ~plantInfo.isUnstable
                score = score + 10;  % Bonus for RHP zeros in stable plants (loop shaping handles well)
            end
            if specialConditions.rhpZeroNearUnstablePole
                score = score - 30;  % Strong penalty for RHP zeros near unstable poles
            end
            if specialConditions.highlyOscillatory && ~plantInfo.isUnstable
                score = score + 5;   % Small bonus for oscillatory modes in stable plants
            end
            
        case 'Pole Placement'
            % Good for state-space models and when precise dynamics are needed
            if isfield(options, 'stateSpace') || isa(G, 'ss')
                score = score + 20;  % Major bonus if state-space available
            else
                score = score - 10;  % Penalty if no state-space model
            end
            
            % Difficulty adjustments
            if strcmp(difficulty, 'Extremely Difficult') && isfield(options, 'stateSpace')
                score = score + 15;  % Bonus for extremely difficult plants if SS available
            elseif strcmp(difficulty, 'Very Difficult') && isfield(options, 'stateSpace')
                score = score + 10;  % Bonus for very difficult plants if SS available
            end
            
            % Special conditions bonuses
            if plantInfo.isUnstable && isfield(options, 'stateSpace')
                score = score + 15;  % Bonus for unstable plants with SS model
            end
            if specialConditions.multipleUnstablePoles && isfield(options, 'stateSpace')
                score = score + 20;  % Extra bonus for multiple unstable poles with SS model
            end
            if specialConditions.highlyOscillatory && isfield(options, 'stateSpace')
                score = score + 15;  % Bonus for oscillatory plants with SS model
            end
            if plantInfo.isHighOrder && isfield(options, 'stateSpace')
                score = score + 10;  % Bonus for high-order plants with SS model
            end
            
        case 'Compensation Controller'
            % Excellent for systems with problematic dynamics
            if strcmp(difficulty, 'Extremely Difficult')
                score = score + 30;  % Major bonus for extremely difficult plants
            elseif strcmp(difficulty, 'Very Difficult')
                score = score + 25;  % Strong bonus for very difficult plants
            elseif strcmp(difficulty, 'Difficult')
                score = score + 20;  % Bonus for difficult plants
            end
            
            % Special conditions bonuses
            if specialConditions.multipleUnstablePoles
                score = score + 25;  % Strong bonus for multiple unstable poles
            end
            if specialConditions.highlyUnstable
                score = score + 20;  % Bonus for highly unstable poles
            end
            if specialConditions.closeRHPZeros
                score = score + 15;  % Bonus for close RHP zeros
            end
            if specialConditions.rhpZeroNearUnstablePole
                score = score + 25;  % Strong bonus for challenging RHP zeros
            end
            if specialConditions.complexRHPZeros
                score = score + 10;  % Bonus for complex RHP zeros
            end
            if plantInfo.isHighOrder
                score = score + 10;  % Bonus for high-order systems
            end
    end
    
    % Controller structure considerations
    switch structure
        case 'P'
            % For P controllers, simple methods might be better
            if ismember(method, {'Ziegler-Nichols (Oscillation)', 'Loop-Shaping'})
                score = score + 5;
            end
        case 'PI'
            % For PI controllers, no special adjustment
        case 'PD'
            % For PD controllers, methods with good derivative handling
            if ismember(method, {'Loop-Shaping', 'Pole Placement'})
                score = score + 5;
            end
        case 'PID'
            % For PID controllers, advanced methods are usually better
            if ismember(method, {'Pole Placement', 'Compensation Controller'})
                score = score + 10;
            elseif ismember(method, {'Loop-Shaping'})
                score = score + 5;
            end
    end
    
    % Goal-specific adjustments
    switch options.goal
        case 'Tracking'
            if ismember(method, {'Pole Placement', 'Loop-Shaping'})
                score = score + 5;
            end
        case 'Disturbance Rejection'
            if ismember(method, {'Compensation Controller', 'Loop-Shaping'})
                score = score + 5;
            end
        case 'Robustness'
            if ismember(method, {'Compensation Controller'})
                score = score + 10;
            elseif ismember(method, {'Loop-Shaping'})
                score = score + 5;
            end
    end
    
    % Bandwidth considerations
    if options.bandwidth > 5
        % For high bandwidth requirements
        if ismember(method, {'Pole Placement', 'Loop-Shaping'})
            score = score + 5;
        end
    elseif options.bandwidth < 0.5
        % For low bandwidth requirements
        if ismember(method, {'Compensation Controller', 'Loop-Shaping'})
            score = score + 5;
        end
    end
    
    % Limit score to sensible range
    score = min(max(score, 0), 100);
end

function options = prepareStateSpaceModel(G, options, plantInfo)
    % PREPARESTATESPACEMODEL Prepare state-space model for advanced controller design methods
    % FIXED: Removed invalid return statement and added proper documentation
    %
    % Inputs:
    %   G         - Plant transfer function
    %   options   - Design options structure
    %   plantInfo - Plant information structure
    %
    % Output:
    %   options - Updated options with state-space model included
    
    try
        % First try standard conversion
        disp('Creating state-space model for advanced controller design...');
        options.stateSpace = ss(G);
        [A, B, C, D] = ssdata(options.stateSpace);
        options.stateMatrices = {A, B, C, D};
        disp('State-space model created successfully');
    catch ME
        disp(['Standard state-space conversion failed: ', ME.message]);
        disp('Attempting enhanced state-space conversion...');
        
        % Try enhanced conversions for difficult plants
        try
            % First attempt balanced realization
            options.stateSpace = getEnhancedStateSpace(G, plantInfo);
            [A, B, C, D] = ssdata(options.stateSpace);
            options.stateMatrices = {A, B, C, D};
            disp('Successfully created enhanced balanced state-space model');
        catch ME2
            disp(['Balanced conversion failed: ', ME2.message]);
            disp('Attempting more robust conversion...');
            
            try
                % Manual conversion from transfer function
                [num, den] = tfdata(G, 'v');
                
                % Convert to companion form
                [A, B, C, D] = tf2ss(num, den);
                
                % Create state-space model
                options.stateSpace = ss(A, B, C, D);
                options.stateMatrices = {A, B, C, D};
                disp('Successfully created manual state-space model');
            catch ME3
                disp(['All conversion attempts failed: ', ME3.message]);
                disp('Will proceed without state-space model');
            end
        end
    end
    
    % No return statement needed; MATLAB automatically returns the
    % modified 'options' structure as defined in the function declaration
end

function [K, details] = designControllerWithMethod(method, G, structure, options, plantInfo)
    % DESIGNCONTROLLERWITHMETHOD Design controller using specified method
    % FIXED: Added plantInfo parameter for consistency and added proper documentation
    %
    % Inputs:
    %   method    - Controller design method
    %   G         - Plant transfer function
    %   structure - Controller structure ('P', 'PI', 'PD', 'PID')
    %   options   - Design options structure
    %   plantInfo - Plant information structure
    %
    % Outputs:
    %   K        - Designed controller transfer function
    %   details  - Design details
    
    switch method
        case 'Ziegler-Nichols (Oscillation)'
            [K, details] = designZieglerNicholsOscillation(G, structure, options, plantInfo);
            
        case 'Ziegler-Nichols (Step)'
            [K, details] = designZieglerNicholsStep(G, structure, options, plantInfo);
            
        case 'Aström'
            [K, details] = designAstrom(G, structure, options, plantInfo);
            
        case 'Loop-Shaping'
            [K, details] = designLoopShaping(G, structure, options, plantInfo);
            
        case 'Pole Placement'
            [K, details] = designPolePlacement(G, structure, options, plantInfo);
            
        case 'Compensation Controller'
            [K, details] = designCompensationController(G, structure, options, plantInfo);
            
        otherwise
            error('Unknown design method: %s', method);
    end
end

function [K, details, success] = attemptFallbackDesign(G, structure, options, plantInfo, plantDifficulty, originalError)
    % ATTEMPTFALLBACKDESIGN Try alternative design methods when primary method fails
    % FIXED: Changed parameter name from 'difficulty' to 'plantDifficulty' for consistency
    %
    % Inputs:
    %   G               - Plant transfer function
    %   structure       - Controller structure ('P', 'PI', 'PD', 'PID')
    %   options         - Design options structure
    %   plantInfo       - Plant information structure
    %   plantDifficulty - Plant difficulty classification
    %   originalError   - Original error message from failed design
    %
    % Outputs:
    %   K        - Fallback controller
    %   details  - Design details
    %   success  - Boolean indicating if fallback was successful
    
    % Initialize success flag
    success = false;
    
    % Create fallback sequence based on plant difficulty and characteristics
    % FIXED: Changed 'difficulty' to 'plantDifficulty'
    switch plantDifficulty
        case {'Extremely Difficult', 'Very Difficult'}
            % For most difficult plants, try these methods in sequence
            fallbackMethods = {'Compensation Controller', 'Pole Placement', 'designEmergencyController'};
        case 'Difficult'
            % For difficult plants
            fallbackMethods = {'Compensation Controller', 'Loop-Shaping', 'Pole Placement'};
        otherwise
            % For standard plants
            fallbackMethods = {'Loop-Shaping', 'Compensation Controller', 'Ziegler-Nichols (Oscillation)'};
    end
    
    % Remove any methods that match the original failed method (to avoid repetition)
    originalMethod = regexp(originalError, 'Primary design method \((.*?)\)', 'tokens');
    if ~isempty(originalMethod) && ~isempty(originalMethod{1})
        originalMethod = originalMethod{1}{1};
        fallbackMethods = fallbackMethods(~strcmp(fallbackMethods, originalMethod));
    end
    
    details = sprintf('Original design method failed: %s\n\nFALLBACK DESIGN ATTEMPT\n', originalError);
    details = [details, '------------------------\n'];
    
    % Try each fallback method in sequence
    for i = 1:length(fallbackMethods)
        try
            disp(['Trying fallback method: ', fallbackMethods{i}]);
            details = [details, sprintf('Attempting fallback method %d: %s\n', i, fallbackMethods{i})];
            
            % Special handling for emergency controller
            if strcmp(fallbackMethods{i}, 'designEmergencyController')
                [K, details_method] = designEmergencyController(G, structure, options, plantInfo);
            else
                [K, details_method] = designControllerWithMethod(fallbackMethods{i}, G, structure, options, plantInfo);
            end
            
            details = [details, details_method];
            
            % Verify closed-loop stability
            try
                T = feedback(G*K, 1);
                cl_poles = pole(T);
                is_stable = all(real(cl_poles) < 0);
                
                if is_stable
                    details = [details, sprintf('\nFallback method %s successful!\n', fallbackMethods{i})];
                    success = true;
                    return;
                else
                    details = [details, sprintf('\nFallback method %s produced unstable controller\n', fallbackMethods{i})];
                    
                    % Try gain reduction for stability
                    [num, den] = tfdata(K, 'v');
                    
                    for scale = [0.5, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001]
                        K_test = tf(num * scale, den);
                        T_test = feedback(G * K_test, 1);
                        
                        if all(real(pole(T_test)) < 0)
                            K = K_test;
                            details = [details, sprintf('System stabilized with gain factor: %.4f\n', scale)];
                            success = true;
                            return;
                        end
                    end
                    
                    details = [details, 'Gain adjustment failed to stabilize system\n'];
                end
            catch ME
                details = [details, sprintf('Error verifying stability: %s\n', ME.message)];
            end
            
        catch ME
            details = [details, sprintf('Fallback method %s failed: %s\n', fallbackMethods{i}, ME.message)];
            disp(['Fallback method failed: ', ME.message]);
        end
    end
    
    details = [details, 'All fallback methods failed. Proceeding to emergency controller.\n'];
    success = false;
end

function [K, details] = createEmergencyController(G_tf, structure, options, plantInfo, plantDifficulty)
    % CREATEEMERGENCYCONTROLLER Create a minimal stabilizing controller as last resort
    % FIXED: Changed parameter name from 'difficulty' to 'plantDifficulty' for consistency
    %
    % Inputs:
    %   G_tf           - Plant transfer function
    %   structure      - Controller structure ('P', 'PI', 'PD', 'PID')
    %   options        - Design options structure
    %   plantInfo      - Plant information structure
    %   plantDifficulty - Plant difficulty classification
    %
    % Outputs:
    %   K        - Emergency controller
    %   details  - Design details
    
    details = 'EMERGENCY CONTROLLER DESIGN\n';
    details = [details, '-------------------------\n'];
    details = [details, 'All standard methods failed. Creating conservative emergency controller.\n'];
    
    if plantInfo.isUnstable
        details = [details, 'Plant is unstable. Using specialized stabilization approach.\n'];
        
        try
            % For unstable plants, try progressive stabilization approaches
            [K, emergency_details] = designEmergencyControllerForUnstable(G_tf, structure, options, plantInfo);
            details = [details, emergency_details];
            
            % Verify emergency controller
            T = feedback(G_tf*K, 1);
            if all(real(pole(T)) < 0)
                details = [details, 'Emergency controller successfully stabilizes the system.\n'];
                return;
            else
                details = [details, 'Emergency controller fails to stabilize. Trying one more approach.\n'];
            end
        catch ME
            details = [details, sprintf('Emergency controller error: %s\n', ME.message)];
        end
    end
    
    % Ultra-conservative controller as last resort
    try
        switch structure
            case 'P'
                K = tf(0.001, 1);
                details = [details, 'Created ultra-conservative P controller.\n'];
            case 'PI'
                K = tf([0.001, 0.0001], [1, 0]);
                details = [details, 'Created ultra-conservative PI controller.\n'];
            case 'PD'
                K = tf([0.005, 0.001], [0.001, 1]);
                details = [details, 'Created ultra-conservative PD controller.\n'];
            case 'PID'
                K = tf([0.005, 0.001, 0.0001], [0.001, 1, 0]);
                details = [details, 'Created ultra-conservative PID controller.\n'];
            otherwise
                K = tf(0.001, 1);
                details = [details, 'Created ultra-conservative P controller.\n'];
        end
        
        % Try to verify stability
        try
            T = feedback(G_tf*K, 1);
            if all(real(pole(T)) < 0)
                details = [details, 'Ultra-conservative controller is stable.\n'];
            else
                details = [details, 'Even ultra-conservative controller is unstable. System may be uncontrollable.\n'];
            end
        catch
            details = [details, 'Could not verify stability of emergency controller.\n'];
        end
    catch ME
        % Absolute last resort
        details = [details, sprintf('Error in ultra-conservative controller: %s\n', ME.message)];
        K = tf(0.0001, 1);
        details = [details, 'Created minimal gain controller as absolute last resort.\n'];
    end
    
    details = [details, 'WARNING: Emergency controller used - performance will be very limited.\n'];
    details = [details, 'Manual controller tuning is strongly recommended for this system.\n'];
end

function [K, details] = postProcessController(K, G, details, structure, options, plantInfo)
    % Perform post-processing on the controller to ensure stability and proper structure
    
    % Add post-processing section to details
    details = [details, '\nPOST-PROCESSING CONTROLLER\n'];
    details = [details, '-------------------------\n'];
    
    % 1. Verify and fix controller stability
    try
        [num, den] = tfdata(K, 'v');
        K_poles = roots(den);
        
        if any(real(K_poles) > 0)
            details = [details, 'Controller has unstable poles. Applying stabilization...\n'];
            
            % Stabilize controller poles by reflecting unstable poles
            for i = 1:length(K_poles)
                if real(K_poles(i)) > 0
                    K_poles(i) = -real(K_poles(i)) + imag(K_poles(i))*1i;
                end
            end
            
            % Create new controller with stabilized poles
            den_stable = poly(K_poles);
            K = tf(num, den_stable);
            details = [details, 'Controller poles have been stabilized.\n'];
        else
            details = [details, 'Controller has stable poles - no stabilization needed.\n'];
        end
    catch ME
        details = [details, sprintf('Controller pole analysis failed: %s\n', ME.message)];
    end
    
    % 2. Verify and fix closed-loop stability
    try
        T = feedback(G*K, 1);
        cl_poles = pole(T);
        
        if any(real(cl_poles) > 0)
            details = [details, 'Closed-loop system is unstable. Attempting gain adjustment...\n'];
            
            % Try reducing gain until stable
            [num, den] = tfdata(K, 'v');
            stabilized = false;
            
            for scale = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001, 0.0001]
                K_test = tf(num * scale, den);
                T_test = feedback(G * K_test, 1);
                
                try
                    cl_poles_test = pole(T_test);
                    
                    if all(real(cl_poles_test) < 0)
                        K = K_test;
                        stabilized = true;
                        details = [details, sprintf('System stabilized with gain factor: %.5f\n', scale)];
                        break;
                    end
                catch
                    % Skip this scale factor if pole computation fails
                    continue;
                end
            end
            
            if ~stabilized
                details = [details, 'WARNING: Could not stabilize system with gain adjustment.\n'];
                details = [details, 'The plant may be uncontrollable with this controller structure.\n'];
            end
        else
            details = [details, 'Closed-loop system is stable - no gain adjustment needed.\n'];
        end
    catch ME
        details = [details, sprintf('Closed-loop stability analysis failed: %s\n', ME.message)];
    end
    
    % 3. Ensure the controller has the requested structure
    try
        % Get controller structure
        actual_structure = determineControllerType(K);
        
        if ~strcmpi(actual_structure, structure)
            details = [details, sprintf('Controller has structure %s instead of requested %s. Adjusting...\n', actual_structure, structure)];
            
            % Convert controller to requested structure
            [K_adjusted, adjusted] = adjustControllerStructure(K, structure, options, plantInfo);
            
            if adjusted
                K = K_adjusted;
                details = [details, 'Successfully adjusted controller to requested structure.\n'];
            else
                details = [details, 'Could not adjust controller structure while maintaining stability.\n'];
                details = [details, 'Using original controller structure for stability reasons.\n'];
            end
        else
            details = [details, sprintf('Controller already has the requested %s structure.\n', structure)];
        end
    catch ME
        details = [details, sprintf('Controller structure analysis failed: %s\n', ME.message)];
    end
    
    % 4. Apply numerical conditioning
    try
        [num, den] = tfdata(K, 'v');
        [num_c, den_c] = conditionControllerCoefficients(num, den);
        
        if ~isequal(num, num_c) || ~isequal(den, den_c)
            K = tf(num_c, den_c);
            details = [details, 'Applied numerical conditioning to controller coefficients.\n'];
        else
            details = [details, 'Controller coefficients are already well-conditioned.\n'];
        end
    catch ME
        details = [details, sprintf('Numerical conditioning failed: %s\n', ME.message)];
    end
    
    % 5. Final controller details
    try
        [num, den] = tfdata(K, 'v');
        details = [details, sprintf('\nFINAL CONTROLLER K(s):\n')];
        details = [details, sprintf('Numerator: [%s]\n', mat2str(num, 5))];
        details = [details, sprintf('Denominator: [%s]\n', mat2str(den, 5))];
        
        % Add controller type
        controllerType = determineControllerType(K);
        details = [details, sprintf('Controller type: %s\n', controllerType)];
        
        % Extract key parameters if standard PID form
        if ismember(controllerType, {'P', 'PI', 'PD', 'PID'})
            [Kp, Ki, Kd, Tf] = extractPIDParameters(K);
            
            details = [details, sprintf('Key parameters:\n')];
            if ~isnan(Kp), details = [details, sprintf('Kp = %.5f\n', Kp)]; end
            if ~isnan(Ki), details = [details, sprintf('Ki = %.5f\n', Ki)]; end
            if ~isnan(Kd), details = [details, sprintf('Kd = %.5f\n', Kd)]; end
            if ~isnan(Tf), details = [details, sprintf('Tf = %.5f (derivative filter time constant)\n', Tf)]; end
        end
    catch ME
        details = [details, sprintf('Final controller analysis failed: %s\n', ME.message)];
    end
end

function [K_adjusted, success] = adjustControllerStructure(K, target_structure, options, plantInfo)
    % Adjust controller to match requested structure
    
    % Initialize success flag
    success = false;
    
    % Get current structure
    current_structure = determineControllerType(K);
    
    % If already matching, return
    if strcmpi(current_structure, target_structure)
        K_adjusted = K;
        success = true;
        return;
    end
    
    % Extract transfer function data
    [num, den] = tfdata(K, 'v');
    
    try
        % Extract approximate PID parameters from current controller
        [Kp, Ki, Kd, Tf] = extractPIDParameters(K);
        
        % Apply defaults for missing parameters
        if isnan(Kp), Kp = 1; end
        if isnan(Ki), Ki = 0.1; end
        if isnan(Kd), Kd = 0.1; end
        if isnan(Tf), Tf = options.epsilon; end
        
        % Create new controller with target structure
        switch target_structure
            case 'P'
                K_adjusted = tf(Kp, 1);
                
            case 'PI'
                % If no integral action exists, use conservative Ki
                if isnan(Ki) || Ki == 0
                    Ki = Kp * 0.1;
                end
                K_adjusted = tf([Kp, Ki], [1, 0]);
                
            case 'PD'
                % If no derivative action exists, use conservative Kd
                if isnan(Kd) || Kd == 0
                    Kd = Kp * 0.1;
                end
                K_adjusted = tf([Kd, Kp], [Tf, 1]);
                
            case 'PID'
                % If parameters are missing, use conservative values
                if isnan(Ki) || Ki == 0
                    Ki = Kp * 0.1;
                end
                if isnan(Kd) || Kd == 0
                    Kd = Kp * 0.1;
                end
                K_adjusted = tf([Kd, Kp, Ki], [Tf, 1, 0]);
                
            otherwise
                % Unknown structure, return original
                K_adjusted = K;
                success = false;
                return;
        end
        
        success = true;
        
    catch ME
        % If adjustment fails, return original controller
        K_adjusted = K;
        success = false;
    end
end

function [num_c, den_c] = conditionControllerCoefficients(num, den)
    % Condition controller coefficients for numerical stability
    
    % Remove very small coefficients (numerical noise)
    num_c = num;
    den_c = den;
    
    % Threshold for small coefficients
    threshold = 1e-10;
    
    % Clean numerator
    max_num = max(abs(num_c));
    num_c(abs(num_c) < threshold * max_num) = 0;
    
    % Clean denominator
    max_den = max(abs(den_c));
    den_c(abs(den_c) < threshold * max_den) = 0;
    
    % Scale coefficients if very large or small
    max_coeff = max(max(abs(num_c)), max(abs(den_c)));
    
    if max_coeff > 1e6
        scale_factor = 1e6 / max_coeff;
        num_c = num_c * scale_factor;
    elseif max_coeff < 1e-6
        scale_factor = 1e-6 / max_coeff;
        num_c = num_c * scale_factor;
    end
end

function [Kp, Ki, Kd, Tf] = extractPIDParameters(K)
    % Extract PID parameters from controller transfer function
    
    % Initialize parameters as NaN
    Kp = NaN;
    Ki = NaN;
    Kd = NaN;
    Tf = NaN;
    
    % Get controller type
    controllerType = determineControllerType(K);
    
    % Extract based on controller type
    [num, den] = tfdata(K, 'v');
    
    switch controllerType
        case 'P'
            if length(num) == 1
                Kp = num(1);
                Ki = 0;
                Kd = 0;
                Tf = 0;
            end
            
        case 'PI'
            if length(num) == 2 && length(den) == 2 && den(2) == 0
                Kp = num(1);
                Ki = num(2);
                Kd = 0;
                Tf = 0;
            end
            
        case 'PD'
            if length(num) == 2 && length(den) == 2 && den(2) > 0
                Kd = num(1);
                Kp = num(2);
                Ki = 0;
                Tf = den(1) / num(1);  % Filter time constant
            end
            
        case 'PID'
            if length(num) == 3 && length(den) == 3 && den(3) == 0
                Kd = num(1);
                Kp = num(2);
                Ki = num(3);
                Tf = den(1) / num(1);  % Filter time constant
            end
    end
end

function controllerType = determineControllerType(K)
    % Determine controller type from transfer function
    
    % Get numerator and denominator
    [num, den] = tfdata(K, 'v');
    
    % Remove leading zeros
    num = num(find(abs(num) > 1e-10, 1):end);
    den = den(find(abs(den) > 1e-10, 1):end);
    
    % Check for different controller types
    has_integrator = any(abs(den) < 1e-10);
    
    % Determine type based on transfer function structure
    if length(num) == 1 && length(den) == 1
        % Simple proportional controller
        controllerType = 'P';
    elseif has_integrator && length(num) <= 2 && length(den) <= 2
        % PI controller
        controllerType = 'PI';
    elseif ~has_integrator && length(num) >= 2 && length(den) >= 2 && length(num) <= 3 && length(den) <= 3
        % PD controller with filter
        controllerType = 'PD';
    elseif has_integrator && length(num) >= 3 && length(den) >= 2 && length(den) <= 3
        % PID controller
        controllerType = 'PID';
    elseif length(num) > 3 || length(den) > 3
        % Higher-order controller
        controllerType = 'Higher-Order';
    else
        % Default to unknown type
        controllerType = 'Custom';
    end
end

function plantInfo = createBasicPlantInfo(G)
    % Create basic plantInfo when analyzePlant fails
    
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
    
    % Add placeholder for other fields
    plantInfo.controllability = struct('rank', NaN, 'size', NaN, 'condition', NaN);
    plantInfo.observability = struct('rank', NaN, 'size', NaN, 'condition', NaN);
    plantInfo.modes = struct('eigenvalues', [], 'controllability', [], 'observability', []);
    plantInfo.poorlyControllableModes = [];
    plantInfo.poorlyObservableModes = [];
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
    
    % Get poles and zeros
    p = plantInfo.poles;
    z = plantInfo.zeros;
    
    % Extract unstable poles
    unstable_poles = p(real(p) > 0);
    
    % Sort by real part (most unstable first)
    [~, idx] = sort(real(unstable_poles), 'descend');
    unstable_poles = unstable_poles(idx);
    
    % Initialize compensator
    num_K = 1;
    den_K = 1;
    
    details = [details, sprintf('Plant has %d unstable pole(s)\n', length(unstable_poles))];
    
    % Direct pole cancellation approach
    details = [details, '\n1. Attempting Direct Pole-Zero Cancellation:\n'];
    
    for i = 1:length(unstable_poles)
        pole_i = unstable_poles(i);
        
        if imag(pole_i) ~= 0
            % Skip conjugate pairs, we'll add both together
            if i < length(unstable_poles) && abs(pole_i - conj(unstable_poles(i+1))) < 1e-6
                continue;
            end
            
            if imag(pole_i) > 0
                % For complex poles, create quadratic terms
                real_part = real(pole_i);
                imag_part = imag(pole_i);
                
                details = [details, sprintf('   - Cancelling complex pole at %.3f+%.3fj\n', real_part, imag_part)];
                
                % Create zeros at unstable poles
                quad_term = [1, -2*real_part, real_part^2 + imag_part^2];
                
                % Create stable poles at reflected locations with extra damping
                stable_quad = [1, -4*real_part, 5*(real_part^2 + imag_part^2)];
                
                num_K = conv(num_K, quad_term);
                den_K = conv(den_K, stable_quad);
            end
        else
            % For real poles
            details = [details, sprintf('   - Cancelling real pole at %.3f\n', pole_i)];
            
            num_K = conv(num_K, [1, -pole_i]);
            den_K = conv(den_K, [1, -4*pole_i]);  % Move pole far into LHP
        end
    end
    
    % Use very conservative gain for highly unstable systems
    gain_factor = 0.001;
    num_K = num_K * gain_factor;
    
    details = [details, sprintf('   - Applied conservative gain factor: %.5f\n', gain_factor)];
    
    % Create pre-stabilizing controller
    K_prestab = tf(num_K, den_K);
    
    % Verify pre-stabilization
    try
        prestab_cl = feedback(G * K_prestab, 1);
        prestab_poles = pole(prestab_cl);
        prestab_stable = all(real(prestab_poles) < 0);
        
        if prestab_stable
            details = [details, '   - Pre-stabilization successful!\n'];
            
            % For critical systems, just use the pre-stabilizing controller
            if length(unstable_poles) > 1 || real(unstable_poles(1)) > 5
                details = [details, '   - Using pre-stabilizing controller directly due to high instability\n'];
                K = K_prestab;
                return;
            end
            
            % Otherwise, proceed to adding structure on top of pre-stabilization
            details = [details, '\n2. Adding Requested Controller Structure:\n'];
            
            % Add the requested structure on top of pre-stabilization
            switch structure
                case 'P'
                    K = K_prestab;
                    details = [details, '   - Using P structure directly\n'];
                    
                case 'PI'
                    % Very small integral gain to minimize upset
                    Ki = gain_factor * 0.01;
                    num_PI = [1, Ki];
                    den_PI = [1, 0];
                    
                    K = series(tf(num_PI, den_PI), K_prestab);
                    details = [details, sprintf('   - Added integral action with Ki = %.6f\n', Ki)];
                    
                case 'PD'
                    % Filtered derivative action
                    Kd = gain_factor * 0.1;
                    Tf = 0.2;
                    num_PD = [Kd, 1];
                    den_PD = [Tf, 1];
                    
                    K = series(tf(num_PD, den_PD), K_prestab);
                    details = [details, sprintf('   - Added derivative action with Kd = %.6f, Tf = %.3f\n', Kd, Tf)];
                    
                case 'PID'
                    % PID with very small gains
                    Ki = gain_factor * 0.01;
                    Kd = gain_factor * 0.1;
                    Tf = 0.2;
                    num_PID = [Kd, 1, Ki];
                    den_PID = [Tf, 1, 0];
                    
                    K = series(tf(num_PID, den_PID), K_prestab);
                    details = [details, sprintf('   - Added PID action with Ki = %.6f, Kd = %.6f, Tf = %.3f\n', Ki, Kd, Tf)];
                    
                otherwise
                    K = K_prestab;
            end
            
            % Verify combined stability
            try
                combined_cl = feedback(G * K, 1);
                combined_poles = pole(combined_cl);
                combined_stable = all(real(combined_poles) < 0);
                
                if combined_stable
                    details = [details, '   - Combined controller is stable!\n'];
                else
                    details = [details, '   - Combined controller is unstable, falling back to pre-stabilizing controller\n'];
                    K = K_prestab;
                end
            catch
                details = [details, '   - Error verifying combined stability, falling back to pre-stabilizing controller\n'];
                K = K_prestab;
            end
            
        else
            details = [details, '   - Pre-stabilization failed, attempting gain adjustments\n'];
            
            % Try with different gain factors
            for scale = [0.1, 0.01, 0.001, 0.0001, 0.00001]
                K_test = tf(num_K * scale / gain_factor, den_K);
                try
                    test_cl = feedback(G * K_test, 1);
                    test_poles = pole(test_cl);
                    test_stable = all(real(test_poles) < 0);
                    
                    if test_stable
                        K = K_test;
                        details = [details, sprintf('   - Stabilized with gain factor: %.6f\n', scale)];
                        return;
                    end
                catch
                    continue;
                end
            end
            
            details = [details, '   - All gain adjustments failed\n'];
            details = [details, '\n3. Attempting Alternative Stabilization Approach:\n'];
            
            % Try state-space based approach if plant order is manageable
            if length(p) <= 10
                try
                    % Convert to state-space
                    sys_ss = getEnhancedStateSpace(G, plantInfo);
                    
                    % Extract matrices
                    A = sys_ss.A;
                    B = sys_ss.B;
                    C = sys_ss.C;
                    n = size(A, 1);
                    
                    % Create desired poles in left half-plane
                    desired_poles = zeros(1, n);
                    
                    for i = 1:n
                        if i <= length(p) && real(p(i)) > 0
                            % For unstable poles, reflect and add margin
                            desired_poles(i) = -abs(real(p(i))) * 2;
                        else
                            % For other poles, place them at -1, -2, -3...
                            desired_poles(i) = -i;
                        end
                    end
                    
                    % Attempt pole placement
                    try
                        K_state = place(A, B, desired_poles);
                        details = [details, '   - State feedback design successful\n'];
                        
                        % Convert to transfer function
                        Ac = A - B*K_state;
                        Bc = B;
                        Cc = -K_state;
                        Dc = 0;
                        
                        K_ss = ss(Ac, Bc, Cc, Dc);
                        K = tf(K_ss);
                        
                        % Verify stability
                        cl_ss = feedback(G*K, 1);
                        if all(real(pole(cl_ss)) < 0)
                            details = [details, '   - State feedback controller stabilizes the system\n'];
                            return;
                        else
                            details = [details, '   - State feedback approach failed to stabilize\n'];
                        end
                    catch ME
                        details = [details, sprintf('   - Pole placement failed: %s\n', ME.message)];
                    end
                catch ME
                    details = [details, sprintf('   - State-space approach failed: %s\n', ME.message)];
                end
            else
                details = [details, '   - Plant order too high for state-space approach\n'];
            end
            
            % Last resort - extremely conservative controller
            details = [details, '\n4. Creating Ultra-Conservative Controller:\n'];
            
            switch structure
                case 'P'
                    K = tf(0.0001, 1);
                case 'PI'
                    K = tf([0.0001, 0.00001], [1, 0]);
                case 'PD'
                    K = tf([0.0005, 0.0001], [0.01, 1]);
                case 'PID'
                    K = tf([0.0005, 0.0001, 0.00001], [0.01, 1, 0]);
                otherwise
                    K = tf(0.0001, 1);
            end
            
            details = [details, '   - Created ultra-conservative controller as last resort\n'];
            
            % Try to verify stability
            try
                last_cl = feedback(G*K, 1);
                if all(real(pole(last_cl)) < 0)
                    details = [details, '   - Ultra-conservative controller is stable\n'];
                else
                    details = [details, '   - Even ultra-conservative controller is unstable\n'];
                    details = [details, '   - System may require manual controller design\n'];
                end
            catch
                details = [details, '   - Could not verify stability of ultra-conservative controller\n'];
            end
        end
    catch ME
        details = [details, sprintf('Error in stabilization process: %s\n', ME.message)];
        
        % Create ultra-conservative controller as absolute fallback
        K = tf(0.00001, 1);
        details = [details, 'Created minimal gain controller as absolute last resort.\n'];
    end
end