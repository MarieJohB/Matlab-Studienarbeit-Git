function coefs_clean = cleanup_tf_coefficients(coefs)
    % Handle empty case
    if isempty(coefs)
        coefs_clean = 0;
        return;
    end
    
    % Get the maximum coefficient magnitude
    max_coef = max(abs(coefs));
    
    % If all zeros, just return
    if max_coef == 0
        coefs_clean = coefs;
        return;
    end
    
    % Define threshold relative to largest coefficient
    threshold = max_coef * 1e-10;
    
    % Zero out coefficients below threshold
    coefs_clean = coefs;
    for i = 1:length(coefs_clean)
        if abs(coefs_clean(i)) < threshold
            coefs_clean(i) = 0;
        end
    end
    
    % Remove trailing zeros
    while length(coefs_clean) > 1 && coefs_clean(end) == 0
        coefs_clean = coefs_clean(1:end-1);
    end
end