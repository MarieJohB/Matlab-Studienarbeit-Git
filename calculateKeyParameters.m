function [a, b_abs, b_rel, c, tr, ts] = calculateKeyParameters(t, y_total)
    % Stationary value
    a = y_total(end);
    
    % Overshoot (absolute value)
    b_abs = max(y_total) - a;
    
    % Overshoot (relative value)
    b_rel = b_abs / a;
    
    % Decay ratio
    c = a - min(y_total);
    
    % Rise time (90% of a)
    tr = t(find(y_total >= a * 0.9, 1));
    
    % Settling time (2% band)
    ts = t(find(abs(y_total - a) <= 0.02 * a, 1, 'last'));
    
    % Stationäre Abweichung
    
end

