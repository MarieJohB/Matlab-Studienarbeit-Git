function angles = calculateAngles(lReal, lImag)
    % CALCULATEANGLES - Calculate cumulative angle changes from critical point (-1,0)
    % Implementation specifically for Nyquist stability assessment
    
    % Critical point coordinates
    criticalPoint = [-1, 0];
    
    % Number of points on the Nyquist curve
    numPoints = length(lReal);
    
    % Initialize array for storing cumulative angles
    angles = zeros(numPoints, 1);
    
    % Calculate initial angle from critical point to first point
    initialAngle = atan2(lImag(1) - criticalPoint(2), lReal(1) - criticalPoint(1));
    
    % Determine which axis (real or imaginary) is closer to starting point
    % This is crucial for correct stability assessment
    distToRealAxis = min([abs(initialAngle), abs(initialAngle - pi), abs(initialAngle + pi)]);
    distToImagAxis = min([abs(initialAngle - pi/2), abs(initialAngle + pi/2)]);
    
    % Start counting delta phi from the closest axis
    if distToRealAxis <= distToImagAxis
        % Start from real axis
        disp('Starting angle calculation from real axis');
    else
        % Start from imaginary axis
        disp('Starting angle calculation from imaginary axis');
    end
    
    % Variables to track angle changes
    prevAngle = initialAngle;
    cumulativeAngle = 0;
    
    % Calculate cumulative angle change for each point
    for i = 1:numPoints
        % Calculate current angle from critical point
        currentAngle = atan2(lImag(i) - criticalPoint(2), lReal(i) - criticalPoint(1));
        
        % Calculate angular change from previous point
        deltaAngle = currentAngle - prevAngle;
        
        % Handle angle wrap-around (crossing from -pi to pi or vice versa)
        if deltaAngle > pi
            deltaAngle = deltaAngle - 2*pi;
        elseif deltaAngle < -pi
            deltaAngle = deltaAngle + 2*pi;
        end
        
        % Update cumulative angle
        cumulativeAngle = cumulativeAngle + deltaAngle;
        
        % Store in output array
        angles(i) = cumulativeAngle;
        
        % Update previous angle for next iteration
        prevAngle = currentAngle;
    end
end