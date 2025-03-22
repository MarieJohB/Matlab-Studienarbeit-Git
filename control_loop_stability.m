function [control_loop_stability] = control_loop_stability(G, K)
    % Apply characteristic equation
    [num_G, den_G] = tfdata(G, 'v');
    [num_K, den_K] = tfdata(K, 'v');
    
    % Numerators and denominators are multiplied and the result is added
    char_eq = conv(den_G, den_K) + conv(num_G, num_K);
    
    % Create a transfer function with the characteristic equation as denominator
    char_eq_tf = tf(1, char_eq);
    
    % Use the Routh-Hurwitz method directly
    control_loop_stability = checkRootsNegativeReal(char_eq_tf);
    
    if (control_loop_stability == true)
        disp('Control loop is stable');
    else 
        disp('Control loop is not stable');
    end
end