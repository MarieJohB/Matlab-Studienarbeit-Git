function isUnstableCancellation = checkUnstablePoleZeroCancellation(Z, N)
    % Get poles and zeroes
    poles = roots(Z);  % Poles from the denominator polynomial
    zeroes = roots(N); % Zeroes from the numerator polynomial
    
    % Check for unstable pole-zero cancellation
    isUnstableCancellation = false;
    
    for i = 1:length(poles)
        for j = 1:length(zeroes)
            if abs(poles(i) - zeroes(j)) < 1e-6 && real(poles(i)) > 0
                isUnstableCancellation = true;
                return;
            end
        end
    end
end