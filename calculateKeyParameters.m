function [a, b_abs, b_rel, c, tr, ts] = calculateKeyParameters(t, y_total)



% debugging checks
disp('debugging in calculateKeyParameters')
disp('t:');
disp(length(t));
disp('y_total:');
disp(length(y_total));

    % Stationary value
    a = y_total(end);
    
    % Overshoot (absolute value)
    b_abs = max(y_total) - a;
    
    % Overshoot (relative value)
    b_rel = b_abs / a;
    
    % Decay ratio
    c = a - min(y_total);
    
    % Debugging: check if 90% of a are reached within length of t 
    threshold_90 = a * 0.9;
    disp(['Debugging: 90% of stationary value: ', num2str(threshold_90)]);
    if any(y_total >= threshold_90)
        tr_index = find(y_total >= threshold_90, 1);
        tr = t(tr_index);
    else
        tr = NaN; % if 90% of a will not be reached
        disp('90 % of a will not be reached within the time t.');
    end
    
    % Settling time (2% band)
    ts = NaN; % initialize ts in case the 2% won't be reached
    tolerance_2_percent = 0.02 * a;
    disp(['Debugging: 2% settling time: ', num2str(tolerance_2_percent)]);
    if any(abs(y_total - a) <= tolerance_2_percent)
        ts_index = find(abs(y_total - a) <= tolerance_2_percent, 1, 'last');
        ts = t(ts_index);
    end
end

