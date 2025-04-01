function plantType = determinePlantType(plantInfo)
    % DETERMINEPLANTTYPE Categorize the plant into a specific type to guide
    % controller design recommendations
    
    % Default type
    plantType = 'Standard';
    
    % Check for specific plant characteristics
    if plantInfo.isUnstable && plantInfo.hasRHPZeros
        plantType = 'Difficult';
    elseif plantInfo.isUnstable
        plantType = 'Unstable';
    elseif plantInfo.hasRHPZeros
        plantType = 'NonMinimumPhase';
    elseif plantInfo.hasIntegrator && plantInfo.hasDelay
        plantType = 'IntegratingWithDelay';
    elseif plantInfo.hasIntegrator
        plantType = 'Integrating';
    elseif plantInfo.hasDelay
        plantType = 'WithDelay';
    elseif plantInfo.isHighOrder
        plantType = 'HighOrder';
    end
    
    % Add oscillatory classification if relevant
    if ~isnan(plantInfo.resonancePeak.magnitude) && plantInfo.resonancePeak.magnitude > 2.0
        plantType = [plantType 'Oscillatory'];
    end
end