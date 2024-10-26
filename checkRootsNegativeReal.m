function is_stable = checkRootsNegativeReal(tf_sys)
    % Extract the denominator coefficients
    [~, den] = tfdata(tf_sys, 'v');
    
    % Compute the roots of the denominator polynomial
    den_roots = roots(den);
    
    % Check if all roots have negative real parts
    is_stable = all(real(den_roots) < 0);
end