function is_stable = routh_hurwitz_tf(G)
    % ROUTH_HURWITZ_TF - Stability analysis using Routh-Hurwitz criterion
    %
    % This function determines if a transfer function is asymptotically stable
    % by applying the Routh-Hurwitz criterion to the denominator polynomial.
    %
    % Input:
    %   G - Transfer function model
    %
    % Output:
    %   is_stable - Boolean value (true if stable, false otherwise)
    
    % Extract denominator coefficients
    [~, den] = tfdata(G, 'v');
    
    % Check for NaN coefficients
    if any(isnan(den))
        disp('System contains NaN coefficients - considered unstable');
        is_stable = false;
        return;
    end
    
    % Remove leading zeros if present
    while length(den) > 1 && den(1) == 0
        den = den(2:end);
    end
    
    % Check if denominator is empty after removing zeros
    if isempty(den)
        disp('Invalid transfer function denominator');
        is_stable = false;
        return;
    end
    
    % Normalize denominator by the first coefficient
    den = den / den(1);
    
    % Get the polynomial order
    polyOrder = length(den) - 1;
    
    % If order is 0, the system is stable
    if polyOrder == 0
        is_stable = true;
        return;
    end
    
    % Calculate number of elements in first and second rows
    evenTerms = length(den(1:2:end));
    oddTerms = length(den(2:2:end));
    
    % Number of columns needed for the table
    numColumns = max(evenTerms, oddTerms);
    
    % Initialize Routh-Hurwitz table with zeros
    rhTable = zeros(polyOrder + 1, numColumns);
    
    % Fill first row with coefficients of even powers of s
    rhTable(1, 1:evenTerms) = den(1:2:end);
    
    % Fill second row with coefficients of odd powers of s
    rhTable(2, 1:oddTerms) = den(2:2:end);
    
    % Small value for handling zero entries
    epsilon = 1e-10;
    
    % Flag for marginal stability detection
    marginally_stable = false;
    
    % Compute the rest of the Routh-Hurwitz table
    for i = 3:polyOrder+1
        % Special case: row of all zeros (indicates pairs of symmetric roots)
        if all(abs(rhTable(i-1, :)) < epsilon)
            marginally_stable = true;
            % Use the auxiliary polynomial approach
            order = polyOrder - i + 2;
            temp = zeros(1, numColumns);
            for j = 1:numColumns
                if order >= 0
                    temp(j) = order * rhTable(i-2, j);
                    order = order - 2;
                end
            end
            rhTable(i-1, :) = temp;
        end
        
        % Check for zero in first column
        if abs(rhTable(i-1, 1)) < epsilon
            marginally_stable = true;
            rhTable(i-1, 1) = epsilon; % Replace with a small value
        end
        
        % Compute elements for current row
        for j = 1:numColumns-1
            rhTable(i, j) = (rhTable(i-1, 1) * rhTable(i-2, j+1) - rhTable(i-2, 1) * rhTable(i-1, j+1)) / rhTable(i-1, 1);
        end
    end
    
    % Count sign changes in the first column (indicating RHP poles)
    signChanges = 0;
    for i = 1:polyOrder
        if rhTable(i, 1) * rhTable(i+1, 1) < 0
            signChanges = signChanges + 1;
        elseif abs(rhTable(i+1, 1)) < epsilon
            marginally_stable = true;
        end
    end
    
    % Display Routh-Hurwitz Table (optional)
    disp('Routh-Hurwitz Table:');
    disp(rhTable);
    
    % Determine stability and display result
    if signChanges == 0 && ~marginally_stable
        disp('System is asymptotically stable');
        is_stable = true;
    elseif signChanges == 0 && marginally_stable
        disp('System is marginally stable (poles on imaginary axis) - considered unstable');
        is_stable = false;
    else
        disp('System is unstable');
        is_stable = false;
    end
    
    fprintf('Number of right half-plane poles: %d\n', signChanges);
    
    if marginally_stable
        disp('Warning: System has poles on the imaginary axis');
    end
end