function plantInfo = analyzePlant(G)
    % ANALYZEPLANT Analyze plant characteristics to guide controller design
    % Provides enhanced support for both transfer function and state-space models
    % 
    % Inputs:
    %   G - Plant model (transfer function or state-space)
    %
    % Outputs:
    %   plantInfo - Structure with comprehensive plant analysis information
    
    plantInfo = struct();
    
    % Check if G is a state-space model
    isStateSpace = isa(G, 'ss');
    
    % Get poles and zeros
    if isStateSpace
        % For state-space models, get poles from eigenvalues of A
        plantInfo.poles = eig(G.A);
        
        % Get transmission zeros (more complex for state-space)
        try
            % Use transmission zeros function for state-space models
            sys_tf = tf(G); % First convert to transfer function
            plantInfo.zeros = tzero(G);
            
            % If zero calculation fails, get zeros from transfer function
            if isempty(plantInfo.zeros)
                plantInfo.zeros = zero(sys_tf);
            end
        catch
            % Default to empty if cannot compute zeros
            plantInfo.zeros = [];
        end
    else
        % For transfer functions, use standard methods
        plantInfo.poles = pole(G);
        try
            plantInfo.zeros = zero(G);
        catch
            plantInfo.zeros = [];
        end
    end
    
    % Check stability
    plantInfo.isUnstable = any(real(plantInfo.poles) > 0);
    
    % Check for integrators (poles at the origin)
    plantInfo.hasIntegrator = any(abs(plantInfo.poles) < 1e-6);
    
    % Check for non-minimum phase zeros (RHP zeros)
    plantInfo.hasRHPZeros = any(real(plantInfo.zeros) > 0);
    
    % Check for time delay
    if isStateSpace
        % State-space models typically don't represent pure delays directly
        % Check for Padé-like approximations in high-order models
        plantInfo.hasDelay = false;
        
        % Check if model has alternating signs in characteristic modes
        % which might indicate a Padé approximation
        if size(G.A, 1) >= 4
            [~, p, ~] = zpkdata(tf(G), 'v');
            if length(p) >= 4
                % Check for alternating sign pattern in consecutive coefficients
                [num, den] = tfdata(tf(G), 'v');
                
                % Safely check for alternating signs in numerator and denominator
                hasAlternatingNum = false;
                hasAlternatingDen = false;
                
                if length(num) > 3
                    numSigns = sign(num(1:end-1) .* num(2:end));
                    hasAlternatingNum = any(numSigns < 0);
                end
                
                if length(den) > 3
                    denSigns = sign(den(1:end-1) .* den(2:end));
                    hasAlternatingDen = any(denSigns < 0);
                end
                
                plantInfo.hasDelay = hasAlternatingNum || hasAlternatingDen;
            end
        end
    else
        % Check transfer function for delay characteristics
        [num, den] = tfdata(G, 'v');
        
        % Safely check for alternating signs
        hasAlternatingNum = false;
        hasAlternatingDen = false;
        
        if length(num) > 3
            numSigns = sign(num(1:end-1) .* num(2:end));
            hasAlternatingNum = any(numSigns < 0);
        end
        
        if length(den) > 3
            denSigns = sign(den(1:end-1) .* den(2:end));
            hasAlternatingDen = any(denSigns < 0);
        end
        
        plantInfo.hasDelay = hasAlternatingNum || hasAlternatingDen;
    end
    
    % Determine plant DC gain
    try
        if isStateSpace
            % For state-space, DC gain is C*inv(A)*B + D for stable systems
            if plantInfo.isUnstable || plantInfo.hasIntegrator
                plantInfo.dcGain = Inf;
            else
                plantInfo.dcGain = dcgain(G);
            end
        else
            plantInfo.dcGain = dcgain(G);
        end
    catch
        % For plants with pure integrators or other special cases
        plantInfo.dcGain = Inf;
    end
    
    % Check if high order (more than 2 states)
    if isStateSpace
        plantInfo.isHighOrder = size(G.A, 1) > 2;
    else
        plantInfo.isHighOrder = length(plantInfo.poles) > 2;
    end
    
    % Get step response characteristics if plant is stable
    if ~plantInfo.isUnstable
        try
            t = linspace(0, 100, 1000);
            [y, t] = step(G, t);
            
            % Convert to standardized format if multidimensional
            if size(y, 2) > 1
                y = y(:,1,1); % Take first output, first input
            end
            
            info = stepinfo(y, t);
            plantInfo.stepInfo = info;
            plantInfo.stepResponse = struct('time', t, 'response', y);
            
            % Compute approximate first-order + delay model
            [L, T, K] = estimateFirstOrderPlusDelay(t, y);
            plantInfo.FOPDT = struct('L', L, 'T', T, 'K', K);
        catch
            plantInfo.stepInfo = [];
            plantInfo.stepResponse = [];
            plantInfo.FOPDT = struct('L', NaN, 'T', NaN, 'K', NaN);
        end
    else
        plantInfo.stepInfo = [];
        plantInfo.stepResponse = [];
        plantInfo.FOPDT = struct('L', NaN, 'T', NaN, 'K', NaN);
    end
    
    % Enhanced analysis for state-space models
    if isStateSpace
        % Controllability and observability analysis
        A = G.A;
        B = G.B;
        C = G.C;
        
        % Use SVD for more reliable rank calculation
        Co = ctrb(A, B);
        Ob = obsv(A, C);
        
        sv_ctrb = svd(Co);
        sv_obsv = svd(Ob);
        
        % Set threshold for rank determination
        rank_tol = max(size(A)) * eps(norm(A));
        
        plantInfo.controllability = struct('rank', sum(sv_ctrb > rank_tol), ...
                                          'size', size(A, 1), ...
                                          'condition', max(sv_ctrb)/min(sv_ctrb));
        
        plantInfo.observability = struct('rank', sum(sv_obsv > rank_tol), ...
                                         'size', size(A, 1), ...
                                         'condition', max(sv_obsv)/min(sv_obsv));
        
        % Mode analysis
        [V, D] = eig(A);
        eigA = diag(D);
        
        % Calculate controllability/observability measures for each mode
        mode_ctrb = zeros(length(eigA), 1);
        mode_obsv = zeros(length(eigA), 1);
        
        for i = 1:length(eigA)
            mode_ctrb(i) = norm(B' * V(:,i));
            mode_obsv(i) = norm(C * V(:,i));
        end
        
        % Store mode analysis
        plantInfo.modes = struct('eigenvalues', eigA, ...
                                'controllability', mode_ctrb, ...
                                'observability', mode_obsv);
        
        % Identify poorly controllable/observable modes
        poor_ctrb_idx = find(mode_ctrb < 0.01 * max(mode_ctrb));
        poor_obsv_idx = find(mode_obsv < 0.01 * max(mode_obsv));
        
        plantInfo.poorlyControllableModes = eigA(poor_ctrb_idx);
        plantInfo.poorlyObservableModes = eigA(poor_obsv_idx);
    else
        % For transfer functions, add placeholders for these fields
        plantInfo.controllability = struct('rank', NaN, 'size', NaN, 'condition', NaN);
        plantInfo.observability = struct('rank', NaN, 'size', NaN, 'condition', NaN);
        plantInfo.modes = struct('eigenvalues', [], 'controllability', [], 'observability', []);
        plantInfo.poorlyControllableModes = [];
        plantInfo.poorlyObservableModes = [];
    end
    
    % For unstable plants, estimate stabilizing feedback gain
    if plantInfo.isUnstable
        plantInfo.stabilizingGain = estimateStabilizingGain(G);
    else
        plantInfo.stabilizingGain = 0;
    end
    
    % Add advanced frequency-domain analysis
    try
        % Get frequency response data
        w = logspace(-3, 3, 200);
        [mag, phase] = bode(G, w);
        mag = squeeze(mag);
        phase = squeeze(phase);
        
        % Calculate key frequency-domain metrics
        
        % 1. Try to find resonant peak
        [mag_peak, idx_peak] = max(mag);
        if ~isempty(idx_peak) && idx_peak > 1 && idx_peak < length(w)
            plantInfo.resonancePeak = struct('magnitude', mag_peak, 'frequency', w(idx_peak));
        else
            plantInfo.resonancePeak = struct('magnitude', NaN, 'frequency', NaN);
        end
        
        % 2. Estimate bandwidth (frequency where magnitude drops to 0.707)
        try
            if ~isnan(plantInfo.dcGain) && ~isinf(plantInfo.dcGain)
                ref_mag = 0.707 * abs(plantInfo.dcGain);
            else
                ref_mag = 0.707 * mag(1);
            end
            
            bw_idx = find(mag < ref_mag, 1);
            if ~isempty(bw_idx) && bw_idx > 1
                plantInfo.bandwidth = w(bw_idx);
            else
                plantInfo.bandwidth = NaN;
            end
        catch
            plantInfo.bandwidth = NaN;
        end
        
        % 3. Calculate optimal crossover frequency for loop-shaping
        if ~isnan(plantInfo.bandwidth)
            % For stable systems, typical crossover at 0.5-1x bandwidth
            if ~plantInfo.isUnstable && ~plantInfo.hasRHPZeros
                plantInfo.recommendedCrossover = 0.7 * plantInfo.bandwidth;
            elseif plantInfo.hasRHPZeros
                % For non-minimum phase, crossover should be lower
                min_rhp_zero = min(abs(real(plantInfo.zeros(real(plantInfo.zeros) > 0))));
                plantInfo.recommendedCrossover = min(0.5 * plantInfo.bandwidth, 0.3 * min_rhp_zero);
            else
                % For unstable systems, crossover typically higher than bandwidth
                plantInfo.recommendedCrossover = 1.5 * plantInfo.bandwidth;
            end
        else
            % If bandwidth couldn't be calculated, use pole information
            if plantInfo.isUnstable
                max_unstable_pole = max(real(plantInfo.poles(real(plantInfo.poles) > 0)));
                plantInfo.recommendedCrossover = 2 * max_unstable_pole;
            else
                dom_pole_idx = find(abs(real(plantInfo.poles)) == min(abs(real(plantInfo.poles))));
                if ~isempty(dom_pole_idx)
                    plantInfo.recommendedCrossover = 2 * abs(real(plantInfo.poles(dom_pole_idx(1))));
                else
                    plantInfo.recommendedCrossover = 1.0; % Default
                end
            end
        end
        
        % 4. Identify difficult control frequencies
        plantInfo.difficultFrequencies = [];
        
        % Look for frequencies where phase crosses -180°
        phase_cross_idx = [];
        for i = 1:length(phase)-1
            if (phase(i) > -180 && phase(i+1) < -180) || (phase(i) < -180 && phase(i+1) > -180)
                phase_cross_idx(end+1) = i;
            end
        end
        
        if ~isempty(phase_cross_idx)
            plantInfo.difficultFrequencies = w(phase_cross_idx);
        end
        
        % For RHP zeros, those frequencies are difficult
        if plantInfo.hasRHPZeros
            rhp_zeros = plantInfo.zeros(real(plantInfo.zeros) > 0);
            for i = 1:length(rhp_zeros)
                if imag(rhp_zeros(i)) == 0
                    % Real RHP zero
                    freq = real(rhp_zeros(i));
                    plantInfo.difficultFrequencies(end+1) = freq;
                end
            end
        end
        
        % Sort and remove duplicates
        plantInfo.difficultFrequencies = unique(plantInfo.difficultFrequencies);
        
    catch
        % If frequency analysis fails, initialize with defaults
        plantInfo.resonancePeak = struct('magnitude', NaN, 'frequency', NaN);
        plantInfo.bandwidth = NaN;
        plantInfo.recommendedCrossover = NaN;
        plantInfo.difficultFrequencies = [];
    end
    
    % Calculate relative degree
    try
        [num, den] = tfdata(G, 'v');
        plantInfo.relativeOrder = length(den) - length(num);
    catch
        plantInfo.relativeOrder = NaN;
    end
    
    % Calculate plant type for controller recommendations
    plantInfo.type = determinePlantType(plantInfo);
    
    % Perform additional analysis for non-minimum phase plants
    if plantInfo.hasRHPZeros
        plantInfo.rhpZeroLimitations = analyzeRHPZeroLimitations(plantInfo);
    else
        plantInfo.rhpZeroLimitations = [];
    end
    
    % Calculate default PID parameters as a reference
    plantInfo.defaultPID = calculateDefaultPID(plantInfo);
end







