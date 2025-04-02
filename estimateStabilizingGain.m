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