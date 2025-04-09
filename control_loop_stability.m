function [control_loop_stability] = control_loop_stability(G, K)
    % CONTROL_LOOP_STABILITY Check stability of control loop using Routh-Hurwitz criterion
    %
    % This function determines if a control loop with plant G and controller K
    % is stable by using the existing routh_hurwitz_tf function on the 
    % characteristic equation.
    %
    % Input:
    %   G - Plant transfer function
    %   K - Controller transfer function
    %
    % Output:
    %   control_loop_stability - Boolean indicating loop stability
    
    % Extract numerators and denominators
    [num_G, den_G] = tfdata(G, 'v');
    [num_K, den_K] = tfdata(K, 'v');
    
    % Calculate characteristic equation: den_G * den_K + num_G * num_K
    % This corresponds to 1 + G*K = 0 for the closed-loop characteristic equation
    char_eq = conv(den_G, den_K) + conv(num_G, num_K);
    
    % Create a transfer function with the characteristic equation as denominator
    % and a constant 1 as numerator
    char_eq_tf = tf(1, char_eq);
    
    % Apply Routh-Hurwitz stability criterion using the existing function
    control_loop_stability = routh_hurwitz_tf(char_eq_tf);
    
    % Display result
    if control_loop_stability
        disp('Control loop is stable');
    else
        disp('Control loop is not stable');
    end
end