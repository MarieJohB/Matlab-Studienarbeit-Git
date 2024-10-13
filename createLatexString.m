function latex_str = createLatexString(tf_sys, name, stability, type)
    % Create LaTeX string for the transfer function
    [num, den] = tfdata(tf_sys, 'v');
    num_str = poly2str(num, 's');
    den_str = poly2str(den, 's');
    stability_str = 'Stable';
    if ~stability
        stability_str = 'Unstable';
    end
    latex_str = ['$' name '(s) = \frac{' num_str '}{' den_str '} \Rightarrow \mathrm{' stability_str '} \, | \, ' type '$'];
end