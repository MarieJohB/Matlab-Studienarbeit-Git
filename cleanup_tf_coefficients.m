function coefs_clean = cleanup_tf_coefficients(coefs)
    % Helper function to clean up very small coefficients that might be numerical noise
    % This improves the quality of the resulting transfer function
    
    % Define threshold relative to the magnitude of largest coefficient
    max_coef = max(abs(coefs));
    threshold = max_coef * 1e-10;
    
    % Zero out small coefficients
    coefs_clean = coefs;
    small_indices = abs(coefs) < threshold;
    coefs_clean(small_indices) = 0;
end