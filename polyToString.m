function str = polyToString(poly)
    % Convert polynomial coefficients to a LaTeX string using symbolic math
    syms s;
    polynomial = poly2sym(poly, s);
    str = latex(polynomial);
end
