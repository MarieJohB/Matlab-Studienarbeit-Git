function plantInfo = analyzePlant(G)
    % ANALYZEPLANT Analyze plant characteristics to guide controller design
    % Provides enhanced support for both transfer function and state-space models
    % 
    % Inputs:
    %   G - Plant model (transfer function or state-space)
    %
    % Outputs:
    %   plantInfo - Structure with plant analysis information
    
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
end

function stabilizingGain = estimateStabilizingGain(G)
    % ESTIMATESTABILIZINGGAIN Estimate a stabilizing feedback gain for unstable plants
    
    isStateSpace = isa(G, 'ss');
    
    if isStateSpace
        % For state-space models, use eigenvalues directly
        p = eig(G.A);
    else
        % For transfer functions, use poles
        p = pole(G);
    end
    
    % Find the most unstable pole
    [maxRealPart, idx] = max(real(p));
    
    if maxRealPart <= 0
        % Plant is stable
        stabilizingGain = 0;
        return;
    end
    
    % For a dominant real unstable pole, we need a proportional gain
    % that moves it to the LHP. Use a simple mapping formula.
    if imag(p(idx)) == 0
        % Pure real pole - simple proportional stabilization
        stabilizingGain = maxRealPart * 1.5;
    else
        % Complex pole - use a more conservative gain
        stabilizingGain = maxRealPart * 2;
    end
    
    % For state-space models, check if B vector has suitable magnitude
    % If B is very small, we need a larger gain
    if isStateSpace
        % Get the eigenvectors corresponding to unstable poles
        [V, D] = eig(G.A);
        unstable_idx = find(real(diag(D)) > 0);
        
        % Calculate input coupling to unstable modes
        input_coupling = zeros(length(unstable_idx), 1);
        for i = 1:length(unstable_idx)
            input_coupling(i) = norm(G.B' * V(:,unstable_idx(i)));
        end
        
        % If input coupling is small, increase the gain accordingly
        if ~isempty(input_coupling) && min(input_coupling) < 0.01
            coupling_factor = 0.01 / min(input_coupling);
            stabilizingGain = stabilizingGain * coupling_factor;
        end
    end
end

function [L, T, K] = estimateFirstOrderPlusDelay(t, y)
    % ESTIMATEFIRSTORDERPLUSDELAY Estimate a FOPDT model from step response
    
    % Get steady-state gain
    y_final = y(end);
    K = y_final;
    
    if abs(K) < 1e-6
        % No steady-state response
        L = 0.1;
        T = 1.0;
        return;
    end
    
    % Normalize the response
    y_norm = y / y_final;
    
    % Find 63.2% response point for time constant
    idx_63 = find(y_norm >= 0.632, 1);
    
    if isempty(idx_63)
        % Use 95% of response time as approximation
        idx_95 = find(y_norm >= 0.95, 1);
        if isempty(idx_95)
            T = t(end) / 3;
        else
            T = t(idx_95) / 3;
        end
    else
        T = t(idx_63);
    end
    
    % Find 10% response for delay estimation
    idx_10 = find(y_norm >= 0.1, 1);
    
    if isempty(idx_10)
        L = T * 0.1; % Default approximation
    else
        L = t(idx_10);
    end
    
    % Ensure reasonable values
    L = max(0.01, min(L, T)); % Bound delay to be less than time constant
    T = max(0.1, T);
end