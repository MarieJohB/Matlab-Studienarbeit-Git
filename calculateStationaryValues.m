function [y_inf_1, e_inf_1, y_inf_2, e_inf_2, e_inf_3, cansysjump] = calculateStationaryValues(T_s, S_s, R_s, D1_s, D2_s, G_s)
    % Define the symbolic variable
    syms s;

    % Initialize output variables and convert to symbolic transfer functions
    [T_num, T_den] = tfdata(T_s, 'v');
    [S_num, S_den] = tfdata(S_s, 'v');
    [R_num, R_den] = tfdata(R_s, 'v');
    [D1_num, D1_den] = tfdata(D1_s, 'v');
    [D2_num, D2_den] = tfdata(D2_s, 'v');
    [G_num, G_den] = tfdata(G_s, 'v');

    T_sys = poly2sym(T_num, s) / poly2sym(T_den, s);
    S_sys = poly2sym(S_num, s) / poly2sym(S_den, s);
    R_sys = poly2sym(R_num, s) / poly2sym(R_den, s);
    D1_sys = poly2sym(D1_num, s) / poly2sym(D1_den, s);
    D2_sys = poly2sym(D2_num, s) / poly2sym(D2_den, s);
    G_sys = poly2sym(G_num, s) / poly2sym(G_den, s);

    % Calculate the limits
    y_inf_1 = limit(s * T_sys * R_sys, s, 0);
    e_inf_1 = limit(s * S_sys * R_sys, s, 0);
    y_inf_2 = limit(s * S_sys * D1_sys, s, 0);
    e_inf_2 = limit(-s * S_sys * D1_sys, s, 0);
    e_inf_3 = limit(s * (S_sys * R_sys - G_sys * S_sys * D2_sys), s, 0);
    cansysjump = limit(s * T_sys * R_sys, s, inf);