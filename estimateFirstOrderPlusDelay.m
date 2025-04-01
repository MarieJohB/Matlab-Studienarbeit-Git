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