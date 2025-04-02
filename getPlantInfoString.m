function infoStr = getPlantInfoString(plantInfo)
    % Initialize output string
    infoStr = '';
    
    % Add stability information
    if plantInfo.isUnstable
        infoStr = [infoStr, 'Unstable, '];
    else
        infoStr = [infoStr, 'Stable, '];
    end
    
    % Add integrator information
    if plantInfo.hasIntegrator
        infoStr = [infoStr, 'Has integrator, '];
    end
    
    % Add RHP zeros information
    if plantInfo.hasRHPZeros
        infoStr = [infoStr, 'Non-minimum phase, '];
    end
    
    % Add delay information
    if plantInfo.hasDelay
        infoStr = [infoStr, 'Has delay, '];
    end
    
    % Add order information
    if plantInfo.isHighOrder
        infoStr = [infoStr, 'High-order, '];
    else
        infoStr = [infoStr, 'Low-order, '];
    end
    
    % Add DC gain if available
    if ~isnan(plantInfo.dcGain) && ~isinf(plantInfo.dcGain)
        infoStr = [infoStr, sprintf('DC gain=%.3g', plantInfo.dcGain)];
    else
        infoStr = [infoStr, 'Infinite DC gain'];
    end
end