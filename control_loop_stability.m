function [control_loop_stability] = control_loop_stability(G, K)
%Function to determine if a controll loop is stable

% input: transfer function G(s) and controller K(s) of the system

% applying characteristic quation:
[num_G,den_G] = tfdata(G,'v');
[num_K,den_K] = tfdata(K,'v');

% numerators and denumerators are muliplied and the result is added
char_eq = conv(den_G, den_K) + conv(num_G, num_K);

% now the characteristic quation is checked for stability
% Ruth-Hurwitz-Matrix is used for this
% the input expects a transferfunction
% numerator does not matter, set to one 
char_eq_tf = tf(1, char_eq);

% calling function to check stability
control_loop_stability = routh_hurwitz_tf(char_eq_tf);


if (control_loop_stability == true)
    disp('Control loop is stable');
else 
    disp('Control loop is not stable');
end



% G(s) & K(s) needs to be proper and stable
% S(s) needs to be stable 

%if(proper_tf(G)== true && proper_tf(K) == true && routh_hurwitz_tf(G) == true && routh_hurwitz_tf(K) == true && routh_hurwitz_tf(S) == true )
%    control_loop_stability = true;
%    disp('Control loop is stable');
%else
%    control_loop_stability = false;
%    disp('Control loop is not stable');
%
%end