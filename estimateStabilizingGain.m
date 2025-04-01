function stabilizingGain = estimateStabilizingGain(G)
    % ESTIMATESTABILIZINGGAIN Estimate a stabilizing feedback gain for unstable plants
    
    p = pole(G);
    
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
end