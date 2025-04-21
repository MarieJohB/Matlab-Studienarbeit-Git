function [K, details, score] = robustControllerDesign(G, method, structure, options)
% ROBUSTCONTROLLERDESIGN Comprehensive controller design with enhanced robustness
% Integrates all improved components for automatic controller design
% with particular focus on challenging control plants
%
% Inputs:
%   G         - Plant model (transfer function or state-space)
%   method    - Design method or 'auto' for automatic selection
%               Options: 'auto', 'Compensation Controller', 'Pole Placement',
%                        'Loop-Shaping', 'Ziegler-Nichols (Oscillation)',
%                        'Ziegler-Nichols (Step)', 'Aström'
%   structure - Controller structure: 'P', 'PI', 'PD', 'PID'
%   options   - Design options (see documentation for details)
%
% Outputs:
%   K        - Designed controller as transfer function
%   details  - Text description of design process
%   score    - Performance score (0-100)

% 1. Input validation
if nargin < 3
    error('At least 3 inputs required: G, method, and structure');
end

if nargin < 4
    options = struct();
end

% 2. Perform plant analysis with enhanced diagnostics
try
    disp('Analyzing plant...');
    plantInfo = analyzePlant(G);
    plantInfoString = getPlantInfoString(plantInfo);
    disp(['Plant Analysis: ', plantInfoString]);
    
    % Analyze RHP zero limitations if present
    if plantInfo.hasRHPZeros
        disp('Plant has RHP zeros - analyzing fundamental limitations...');
        rhpLimitations = analyzeRHPZeroLimitations(plantInfo);
        disp(rhpLimitations.description);
    end
    
    % Determine plant difficulty
    [plantDifficulty, specialConditions] = classifyPlantDifficulty(plantInfo);
    disp(['Plant difficulty classification: ', plantDifficulty]);
    
    % Display special warnings for challenging plants
    if strcmp(plantDifficulty, 'Extremely Difficult') || strcmp(plantDifficulty, 'Very Difficult')
        disp('⚠️ WARNING: This plant has severe control challenges');
        disp('   - Achievable performance may be limited');
        disp('   - Consider redesigning the plant if possible');
    end
    
catch ME
    warning('Plant analysis error: %s. Creating basic plant info.', ME.message);
    plantInfo = createBasicPlantInfo(G);
    plantDifficulty = 'Unknown';
    specialConditions = struct();
end

% 3. Apply numerical conditioning to input plant if needed
try
    [num, den] = tfdata(G, 'v');
    [num_c, den_c] = conditionTransferFunction(num, den);
    
    % Check if conditioning was significant
    if ~isequal(num, num_c) || ~isequal(den, den_c)
        disp('Applied numerical conditioning to plant model for better accuracy');
        G_c = tf(num_c, den_c);
    else
        G_c = G;
    end
catch
    G_c = G;
end

% 4. Adjust options based on plant analysis
options = adjustOptionsForDifficulty(options, plantInfo, plantDifficulty, specialConditions);

% 5. Design controller with enhanced error handling
try
    disp(['Designing controller using ', method, ' method with ', structure, ' structure...']);
    

    [K, details, score] = design_controller_auto(G_c, method, structure, options);

    
    % If controller design succeeded but score is very poor, try fallback approach
    if score < 20 && ~isnan(score)
        disp('Controller design succeeded but with very poor score. Trying fallback approach...');
        
        if strcmp(plantDifficulty, 'Extremely Difficult') || strcmp(plantDifficulty, 'Very Difficult')
            % For very difficult plants, use specialized emergency controller
            [K_emergency, details_emergency] = designEmergencyController(G_c, structure, options, plantInfo);
            
            % Test if emergency controller is better
            try
                score_emergency = evaluateController(K_emergency, G_c, options.goal, options.phaseMargin, options.overshoot, options.settlingTime, options.bandwidth, plantInfo);
                
                if isnan(score) || score_emergency > score
                    K = K_emergency;
                    details = ['EMERGENCY CONTROLLER DESIGN (FALLBACK)\n', details_emergency];
                    score = score_emergency;
                    disp('Using emergency controller as it achieved better performance');
                end
            catch
                % Keep original controller if evaluation fails
            end
        else
            % For less difficult plants, try compensation controller
            try
                [K_comp, details_comp] = designCompensationController(G_c, structure, options, plantInfo);
                
                % Test if compensation controller is better
                score_comp = evaluateController(K_comp, G_c, options.goal, options.phaseMargin, options.overshoot, options.settlingTime, options.bandwidth, plantInfo);
                
                if isnan(score) || score_comp > score
                    K = K_comp;
                    details = ['COMPENSATION CONTROLLER DESIGN (FALLBACK)\n', details_comp];
                    score = score_comp;
                    disp('Using compensation controller as it achieved better performance');
                end
            catch
                % Keep original controller if compensation design fails
            end
        end
    end
    
catch ME
    warning('Controller design failed: %s', ME.message);
    disp('Attempting emergency controller design...');
    
    % Emergency fallback for design failure
    try
        [K, details] = designEmergencyController(G_c, structure, options, plantInfo);
        details = ['EMERGENCY CONTROLLER DESIGN (AFTER ERROR)\n', details];
        
        % Evaluate the emergency controller
        try
            score = evaluateController(K, G_c, options.goal, options.phaseMargin, options.overshoot, options.settlingTime, options.bandwidth, plantInfo);
        catch
            score = NaN;
        end
    catch ME2
        error('All controller design methods failed: %s', ME2.message);
    end
end

% 6. Apply numerical robustness improvements
try
    disp('Applying numerical robustness improvements...');
    [K, details] = makeControllerRobust(K, details);
catch ME
    warning('Robustness improvements failed: %s', ME.message);
end

% 7. Final validation of closed-loop stability
try
    disp('Performing final stability check...');
    cl = feedback(G_c*K, 1);
    cl_poles = pole(cl);
    
    if any(real(cl_poles) >= 0)
        warning('Final controller does not stabilize the plant. Attempting gain adjustment...');
        
        % Try scaling the gain to achieve stability
        [num, den] = tfdata(K, 'v');
        stabilized = false;
        
        for scale = [0.5, 0.2, 0.1, 0.05, 0.01, 0.001]
            K_test = tf(num * scale, den);
            cl_test = feedback(G_c * K_test, 1);
            
            if all(real(pole(cl_test)) < 0)
                K = K_test;
                details = [details, sprintf('\n\nFINAL ADJUSTMENT: Gain reduced by factor %.3f for stability', scale)];
                disp(['Stability achieved with gain factor: ', num2str(scale)]);
                stabilized = true;
                
                % Update score with new controller
                try
                    score = evaluateController(K, G_c, options.goal, options.phaseMargin, options.overshoot, options.settlingTime, options.bandwidth, plantInfo);
                catch
                    % Keep previous score if evaluation fails
                end
                
                break;
            end
        end
        
        if ~stabilized
            warning('Unable to stabilize plant with gain adjustment. Plant may be difficult to control.');
            details = [details, '\n\nWARNING: Final controller does not stabilize the plant.'];
            details = [details, '\nConsider redesigning the plant or using advanced control techniques.'];
        end
    else
        disp('Final controller successfully stabilizes the plant.');
        
        % Calculate stability margins for information
        try
            [Gm, Pm, Wcg, Wcp] = margin(G_c*K);
            details = [details, sprintf('\n\nFINAL STABILITY MARGINS:')];
            details = [details, sprintf('\n- Phase Margin: %.2f degrees at %.3f rad/s', Pm, Wcp)];
            details = [details, sprintf('\n- Gain Margin: %.2f dB at %.3f rad/s', 20*log10(Gm), Wcg)];
        catch
            % Skip margin calculation if it fails
        end
    end
catch ME
    warning('Final stability check failed: %s', ME.message);
end

% 8. Add implementation recommendations
details = [details, '\n\nIMPLEMENTATION NOTES:'];

% Add sampling rate recommendation
try
    p = pole(K);
    fastest_pole = max(abs(p));
    recommended_sampling = fastest_pole * 20;
    
    details = [details, sprintf('\n- Recommended minimum sampling rate: %.1f rad/s', recommended_sampling)];
    details = [details, sprintf('\n- Recommended minimum sampling time: %.5f s', 2*pi/recommended_sampling)];
catch
    % Skip if pole calculation fails
end

% Add anti-windup recommendation for controllers with integral action
if strcmpi(structure, 'PI') || strcmpi(structure, 'PID')
    details = [details, '\n- Implement anti-windup protection for integral term'];
    
    % Extract controller parameters for anti-windup recommendation
    try
        [Kp, Ki, Kd, Tf] = extractPIDParameters(K);
        
        if ~isnan(Ki) && ~isnan(Kp) && Ki > 0
            windup_time = 10 * Kp / Ki;
            details = [details, sprintf('\n- Anti-windup tracking time constant: %.2f s', windup_time)];
        end
    catch
        % Skip if parameter extraction fails
    end
end

% 9. Summarize design results
disp('Controller design complete.');
disp(['Controller performance score: ', num2str(score), '/100']);

% Extract key controller parameters for summary
try
    [Kp, Ki, Kd, Tf] = extractPIDParameters(K);
    
    disp('Controller parameters:');
    if ~isnan(Kp), disp(['- Kp = ', num2str(Kp)]); end
    if ~isnan(Ki), disp(['- Ki = ', num2str(Ki)]); end
    if ~isnan(Kd), disp(['- Kd = ', num2str(Kd)]); end
    if ~isnan(Tf), disp(['- Tf = ', num2str(Tf), ' (derivative filter time constant)']); end
catch
    % Skip parameter display if extraction fails
end

return;

end

%% Helper Functions

function options = adjustOptionsForDifficulty(options, plantInfo, difficulty, specialConditions)
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
    
    % Adjust options based on difficulty
    switch difficulty
        case 'Extremely Difficult'
            % Conservative settings for extremely difficult plants
            options.robustness = 'Very High';
            options.damping = min(1.5, options.damping * 1.2);
            options.phaseMargin = max(60, options.phaseMargin);
            options.epsilon = max(0.2, options.epsilon);
            
            % Adjust bandwidth if not explicitly set by user
            if ~options.userSetBandwidth
                if specialConditions.multipleUnstablePoles
                    % Get limit for unstable poles
                    max_unstable_pole = max(real(plantInfo.poles(real(plantInfo.poles) > 0)));
                    options.bandwidth = min(options.bandwidth, max_unstable_pole * 0.2);
                elseif plantInfo.hasRHPZeros
                    % Get limit for RHP zeros
                    min_rhp_zero = min(real(plantInfo.zeros(real(plantInfo.zeros) > 0)));
                    options.bandwidth = min(options.bandwidth, min_rhp_zero * 0.3);
                else
                    options.bandwidth = min(options.bandwidth, 0.3);
                end
                
                disp(['Adjusted bandwidth to ', num2str(options.bandwidth), ' rad/s for extremely difficult plant']);
            end
            
        case 'Very Difficult'
            % More conservative settings for very difficult plants
            options.robustness = 'High';
            options.damping = min(1.2, options.damping * 1.1);
            options.phaseMargin = max(55, options.phaseMargin);
            options.epsilon = max(0.15, options.epsilon);
            
            % Adjust bandwidth if not explicitly set by user
            if ~options.userSetBandwidth
                if specialConditions.multipleUnstablePoles || specialConditions.highlyUnstable
                    max_unstable_pole = max(real(plantInfo.poles(real(plantInfo.poles) > 0)));
                    options.bandwidth = min(options.bandwidth, max_unstable_pole * 0.3);
                elseif plantInfo.hasRHPZeros
                    min_rhp_zero = min(real(plantInfo.zeros(real(plantInfo.zeros) > 0)));
                    options.bandwidth = min(options.bandwidth, min_rhp_zero * 0.4);
                else
                    options.bandwidth = min(options.bandwidth, 0.5);
                end
                
                disp(['Adjusted bandwidth to ', num2str(options.bandwidth), ' rad/s for very difficult plant']);
            end
            
        case 'Difficult'
            % Moderately conservative settings for difficult plants
            options.robustness = 'High';
            options.phaseMargin = max(50, options.phaseMargin);
            options.epsilon = max(0.1, options.epsilon);
            
            % Adjust bandwidth if not explicitly set by user
            if ~options.userSetBandwidth && plantInfo.hasRHPZeros
                min_rhp_zero = min(real(plantInfo.zeros(real(plantInfo.zeros) > 0)));
                options.bandwidth = min(options.bandwidth, min_rhp_zero * 0.5);
                disp(['Adjusted bandwidth to ', num2str(options.bandwidth), ' rad/s due to RHP zero limitations']);
            elseif ~options.userSetBandwidth && plantInfo.isUnstable
                max_unstable_pole = max(real(plantInfo.poles(real(plantInfo.poles) > 0)));
                options.bandwidth = min(options.bandwidth, max_unstable_pole * 0.5);
                disp(['Adjusted bandwidth to ', num2str(options.bandwidth), ' rad/s for unstable plant']);
            end
    end
    
    return;
end

function [Kp, Ki, Kd, Tf] = extractPIDParameters(K)
    % Extract PID parameters from controller transfer function
    
    % Initialize parameters as NaN
    Kp = NaN;
    Ki = NaN;
    Kd = NaN;
    Tf = NaN;
    
    % Get controller type
    [num, den] = tfdata(K, 'v');
    
    % Remove leading zeros
    num = num(find(abs(num) > 1e-10, 1):end);
    den = den(find(abs(den) > 1e-10, 1):end);
    
    % Check for standard forms
    if length(num) == 1 && length(den) == 1
        % P controller
        Kp = num(1) / den(1);
        Ki = 0;
        Kd = 0;
        Tf = 0;
    elseif length(num) == 2 && length(den) == 2 && abs(den(2)) < 1e-10
        % PI controller
        Kp = num(1) / den(1);
        Ki = num(2) / den(1);
        Kd = 0;
        Tf = 0;
    elseif length(num) == 2 && length(den) == 2 && den(2) > 0
        % PD controller with filter
        Kd = num(1) / den(1);
        Kp = num(2) / den(1);
        Ki = 0;
        Tf = den(2) / den(1) / Kd;
    elseif length(num) == 3 && length(den) == 3 && abs(den(3)) < 1e-10
        % PID controller with filter
        Kd = num(1) / den(1);
        Kp = num(2) / den(1);
        Ki = num(3) / den(1);
        Tf = den(2) / den(1) / Kd;
    else
        % Not a standard PID form
        % Try to approximate with frequency response
        try
            w = logspace(-3, 3, 10);
            [mag, phase] = bode(K, w);
            mag = squeeze(mag);
            phase = squeeze(phase);
            
            % Approximate PID parameters
            if phase(1) <= -85
                % Has integral action
                Ki = mag(1) * w(1);  % Approximate Ki from low-frequency response
                
                % Approximate Kp from mid-frequency response
                mid_idx = ceil(length(w)/2);
                Kp = mag(mid_idx);
                
                % Approximate Kd from high-frequency response if derivative action present
                if phase(end) >= 45
                    Kd = mag(end) / w(end);
                else
                    Kd = 0;
                end
            else
                % No integral action
                Ki = 0;
                
                % Approximate Kp from low-frequency response
                Kp = mag(1);
                
                % Approximate Kd from high-frequency response if derivative action present
                if phase(end) >= 45
                    Kd = mag(end) / w(end);
                else
                    Kd = 0;
                end
            end
            
            % Approximate filter constant
            if Kd > 0
                high_idx = length(w);
                expected_mag = Kd * w(high_idx);
                if mag(high_idx) < expected_mag * 0.7
                    Tf = 1 / (w(high_idx) * sqrt(expected_mag / mag(high_idx) - 1));
                else
                    Tf = 0.1;  % Default if cannot determine
                end
            else
                Tf = 0;
            end
        catch
            % Keep as NaN if approximation fails
        end
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
    
    % Add empty fields for consistency
    plantInfo.controllability = struct('rank', NaN, 'size', NaN, 'condition', NaN);
    plantInfo.observability = struct('rank', NaN, 'size', NaN, 'condition', NaN);
    plantInfo.modes = struct('eigenvalues', [], 'controllability', [], 'observability', []);
    plantInfo.poorlyControllableModes = [];
    plantInfo.poorlyObservableModes = [];
end

function info_str = getPlantInfoString(plantInfo)
    % Create a formatted string with plant information
    
    % Initialize output string
    info_str = '';
    
    % Add stability information
    if plantInfo.isUnstable
        info_str = [info_str, 'Unstable, '];
    else
        info_str = [info_str, 'Stable, '];
    end
    
    % Add integrator information
    if plantInfo.hasIntegrator
        info_str = [info_str, 'Has integrator, '];
    end
    
    % Add RHP zeros information
    if plantInfo.hasRHPZeros
        info_str = [info_str, 'Non-minimum phase, '];
    end
    
    % Add delay information
    if plantInfo.hasDelay
        info_str = [info_str, 'Has delay, '];
    end
    
    % Add order information
    if plantInfo.isHighOrder
        info_str = [info_str, 'High-order, '];
    else
        info_str = [info_str, 'Low-order, '];
    end
    
    % Add DC gain if available
    if ~isnan(plantInfo.dcGain) && ~isinf(plantInfo.dcGain)
        info_str = [info_str, sprintf('DC gain=%.3g', plantInfo.dcGain)];
    else
        info_str = [info_str, 'Infinite DC gain'];
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
        min_num = min(abs(num(abs(num) > eps)));
        min_den = min(abs(den(abs(den) > eps)));
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