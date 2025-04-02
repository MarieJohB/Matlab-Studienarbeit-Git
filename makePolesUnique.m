function uniquePoles = makePolesUnique(poles, epsilon)
    % MAKEPOLESUNIQUE - Ensure poles are unique by adding small perturbations
    %
    % This function takes a set of poles and ensures they are uniquely 
    % different from each other by adding small perturbations. This is 
    % necessary for SISO systems where repeated poles can cause issues.
    %
    % Parameters:
    %   poles   - Array of pole locations
    %   epsilon - Optional scale factor for perturbations (default: 1e-6)
    %
    % Returns:
    %   uniquePoles - Array of unique pole locations
    
    if nargin < 2
        epsilon = 1e-6; % Default small perturbation
    end
    
    uniquePoles = poles;
    n = length(poles);
    
    % Check if poles are already unique
    if length(unique(poles)) == length(poles)
        return;
    end
    
    % Group poles by complex conjugate pairs and preserve them
    complexConjugatePairs = [];
    for i = 1:n
        for j = i+1:n
            if abs(real(poles(i)) - real(poles(j))) < epsilon && ...
               abs(imag(poles(i)) + imag(poles(j))) < epsilon && ...
               imag(poles(i)) ~= 0
                % This is a complex conjugate pair - mark for preservation
                complexConjugatePairs(end+1,:) = [i j];
                break; % Once found a pair for pole i, move to next pole
            end
        end
    end
    
    % Identify duplicate poles (excluding complex conjugate pairs)
    duplicateGroups = {};
    alreadyMatched = false(1, n);
    
    for i = 1:n
        if alreadyMatched(i)
            continue;
        end
        
        currentGroup = i;
        
        for j = i+1:n
            if alreadyMatched(j)
                continue;
            end
            
            % Skip checking complex conjugate pairs
            isPart = false;
            for k = 1:size(complexConjugatePairs, 1)
                if (complexConjugatePairs(k, 1) == i && complexConjugatePairs(k, 2) == j) || ...
                   (complexConjugatePairs(k, 1) == j && complexConjugatePairs(k, 2) == i)
                    isPart = true;
                    break;
                end
            end
            
            if isPart
                continue;
            end
            
            % Check for close/duplicate poles
            if abs(poles(i) - poles(j)) < epsilon
                currentGroup = [currentGroup j];
                alreadyMatched(j) = true;
            end
        end
        
        if length(currentGroup) > 1
            duplicateGroups{end+1} = currentGroup;
        end
        
        alreadyMatched(i) = true;
    end
    
    % Process each group of duplicates
    for i = 1:length(duplicateGroups)
        group = duplicateGroups{i};
        for j = 1:length(group)
            idx = group(j);
            % Add perturbation based on position in group
            offset = (j-1) * 2 * epsilon;
            uniquePoles(idx) = uniquePoles(idx) * (1 + offset);
        end
    end
    
    % Restore complex conjugate pairs
    for i = 1:size(complexConjugatePairs, 1)
        idx1 = complexConjugatePairs(i, 1);
        idx2 = complexConjugatePairs(i, 2);
        
        % Ensure they remain conjugates
        realPart = real(uniquePoles(idx1)); 
        imagPart = abs(imag(uniquePoles(idx1)));
        
        uniquePoles(idx1) = realPart + 1i * imagPart;
        uniquePoles(idx2) = realPart - 1i * imagPart;
    end
    
    % Final check to make sure all poles are now unique
    if length(unique(uniquePoles)) < length(uniquePoles)
        % Force uniqueness by adding increasingly larger perturbations
        for i = 1:n
            uniquePoles(i) = uniquePoles(i) * (1 + i * epsilon * 10);
        end
    end
end