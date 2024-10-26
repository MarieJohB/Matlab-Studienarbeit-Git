function is_stable = routh_hurwitz_tf(G)
    % Extract numerator and denominator coefficients
    [num, den] = tfdata(G, 'v');
    coeffVector = den;
    ceoffLength = length(coeffVector);
    rhTableColumn = round(ceoffLength/2);

    % Initialize Routh-Hurwitz table with zeros
    rhTable = zeros(ceoffLength, rhTableColumn);

    % Compute first row of the table
    rhTable(1, :) = coeffVector(1:2:ceoffLength);
    
    % Check if length of coefficients vector is even or odd
    if (rem(ceoffLength, 2) ~= 0)
        % If odd, second row of table will be
        rhTable(2, 1:rhTableColumn - 1) = coeffVector(2:2:ceoffLength);
    else
        % If even, second row of table will be
        rhTable(2, :) = coeffVector(2:2:ceoffLength);
    end
    
    %% Calculate Routh-Hurwitz table's rows
    % Set epss as a small value
    epss = 0.01;
    
    % Calculate other elements of the table
    for i = 3:ceoffLength
        % Special case: row of all zeros
        if all(rhTable(i-1, :) == 0)
            order = ceoffLength - i;
            cnt1 = 0;
            cnt2 = 1;
            for j = 1:rhTableColumn - 1
                rhTable(i-1, j) = (order - cnt1) * rhTable(i-2, cnt2);
                cnt2 = cnt2 + 1;
                cnt1 = cnt1 + 2;
            end
        end
        
        for j = 1:rhTableColumn - 1
            % First element of upper row
            firstElemUpperRow = rhTable(i-1, 1);
            
            % Compute each element of the table
            rhTable(i, j) = ((rhTable(i-1, 1) * rhTable(i-2, j+1)) - ...
                (rhTable(i-2, 1) * rhTable(i-1, j+1))) / firstElemUpperRow;
        end
        
        % Special case: zero in the first column
        if rhTable(i, 1) == 0
            rhTable(i, 1) = epss;
        end
    end
    
    %% Compute number of right hand side poles (unstable poles)
    % Initialize unstable poles with zero
    unstablePoles = 0;
    
    % Check change in signs
    for i = 1:ceoffLength - 1
        if sign(rhTable(i, 1)) * sign(rhTable(i+1, 1)) == -1
            unstablePoles = unstablePoles + 1;
        end
    end
    
    % Display Routh-Hurwitz Table
    disp('Routh-Hurwitz Table:');
    disp(rhTable);
    
    % Display the stability result
    if unstablePoles == 0
        disp('~~~~~> It is a stable system! <~~~~~');
        is_stable = true;
    else
        disp('~~~~~> It is an unstable system! <~~~~~');
        is_stable = false;
    end
    fprintf('\n Number of right hand side poles = %2.0f\n', unstablePoles);
    
end