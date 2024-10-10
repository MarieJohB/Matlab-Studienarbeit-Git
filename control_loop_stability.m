function [control_loop_stability] = control_loop_stability(G,K,S)
%Function to determin if a controll loop is stable 
% G(s) & K(s) needs to be proper and stable
% S(s) needs to be stable 

if(proper_tf(G)== true && proper_tf(K) == true && routh_hurwitz_tf(G) == true && routh_hurwitz_tf(K) == true && routh_hurwitz_tf(S) == true )
    control_loop_stability = true;
    disp('Control loop is stable');
else
    control_loop_stability = false;
    disp('Control loop is not stable');

end