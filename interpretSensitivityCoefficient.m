function interpretation = interpretSensitivityCoefficient(coeff, param, metric)
    % INTERPRETSENSITIVITYCOEFFICIENT - Provide interpretation of sensitivity coefficient
    
    absCoeff = abs(coeff);
    
    % Get appropriate unit for the metric
    switch metric
        case 'Settling Time'
            unit = 'seconds';
        case 'Overshoot'
            unit = '%';
        case 'Rise Time'
            unit = 'seconds';
        case 'Bandwidth'
            unit = 'rad/s';
        case 'Phase Margin'
            unit = 'degrees';
        case 'Gain Margin'
            unit = 'dB';
        otherwise
            unit = 'units';
    end
    
    % Direction of effect
    if coeff > 0
        direction = 'increases';
    elseif coeff < 0
        direction = 'decreases';
    else
        direction = 'does not affect';
    end
    
    % Magnitude of effect
    if absCoeff < 0.1
        strength = 'very little';
    elseif absCoeff < 0.5
        strength = 'slightly';
    elseif absCoeff < 1
        strength = 'moderately';
    elseif absCoeff < 2
        strength = 'significantly';
    else
        strength = 'dramatically';
    end
    
    % Compose the interpretation
    interpretation = sprintf('A 1%% increase in %s %s %s the %s by approximately %.2f%% (%s).', ...
        param, strength, direction, metric, coeff, unit);
    
    % Special case for Gain Margin
    if strcmp(metric, 'Gain Margin')
        % Add specific advice for gain margin
        if coeff < 0
            if absCoeff > 0.5
                interpretation = [interpretation, ' This significant negative sensitivity may pose stability concerns as the parameter increases.'];
            else
                interpretation = [interpretation, ' The negative sensitivity indicates reduced stability margins as the parameter increases.'];
            end
        else
            interpretation = [interpretation, ' The positive sensitivity indicates improved stability margins as the parameter increases.'];
        end
    else
        % Additional insights for other metrics
        if absCoeff > 2
            interpretation = [interpretation, ' The system is highly sensitive to this parameter.'];
        elseif absCoeff < 0.1
            interpretation = [interpretation, ' The system is robust to changes in this parameter.'];
        end
    end
    
    % Add recommendations based on the parameter and coefficient
    if contains(param, 'Gain') && absCoeff > 1
        interpretation = [interpretation, ' Consider detuning this gain for more robust performance.'];
    elseif contains(param, 'Time Constant') && absCoeff > 1
        interpretation = [interpretation, ' The system is highly sensitive to plant time constants. Consider robust control design methods.'];
    end
end