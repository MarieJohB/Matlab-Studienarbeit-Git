function y_t = tranfertosymbolic(sys, t)
    % Convert transfer function to symbolic representation and evaluate it over t
    [num, den] = tfdata(sys, 'v');
    syms s t_sym;
    G_s = poly2sym(num, s) / poly2sym(den, s);
    G_t = ilaplace(G_s, s, t_sym);
    y_t = double(subs(G_t, t_sym, t));
end