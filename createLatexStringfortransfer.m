function latex_str = createLatexStringfortransfer(tf_sys, name)
    % Create LaTeX string for the transfer function
    [num, den] = tfdata(tf_sys, 'v');
    
    % Convert numerator and denominator to strings
    num_str = poly2str(num, 's');
    den_str = poly2str(den, 's');
    
    % Convert scientific notation to LaTeX format
    num_str = regexprep(num_str, 'e([-+]?\d+)', '\\cdot10^{\$1}');
    den_str = regexprep(den_str, 'e([-+]?\d+)', '\\cdot10^{\$1}');
    
    % Create the LaTeX string
    latex_str = ['$' name '(s) = \frac{' num_str '}{' den_str '}$'];
end
