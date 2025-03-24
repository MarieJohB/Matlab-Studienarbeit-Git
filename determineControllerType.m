function controllerType = determineControllerType(K)
    % DETERMINECONTROLLERTYPE - Determine controller type from transfer function
    
    % Get numerator and denominator
    [num, den] = tfdata(K, 'v');
    
    % Remove leading zeros
    num = num(find(num ~= 0, 1):end);
    den = den(find(den ~= 0, 1):end);
    
    % Check for different controller types
    has_integrator = any(abs(den) < 1e-10);
    
    % Determine type based on transfer function structure
    if length(num) == 1 && length(den) == 1
        % Simple proportional controller
        controllerType = 'P';
    elseif has_integrator && length(num) <= 2
        % PI controller
        controllerType = 'PI';
    elseif ~has_integrator && length(num) >= 2 && length(den) >= 2
        % PD controller with filter
        controllerType = 'PD';
    elseif has_integrator && length(num) >= 2
        % PID controller
        controllerType = 'PID';
    else
        % Default to unknown type
        controllerType = 'unknown';
    end
end