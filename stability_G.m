function [ev_G, stab_G, roots_G] = stability_G(G)
% caluculation of eigenvalues
% checking for stability of the controlled systems based on calculation of
% eigenvalues / poles  N_G = 0

% input: G (controlled system)
% output: 
% ev_G: eigenvalues of G
% stab_G: status whether G is stable or not
% roots_G: roots of G


syms s; % set 's' as sybolic variable

% calculation of eigenvalues / poles of G
ev_G = eig(G);

% calculation of poles and roots (Nullstellen) of G
 
[n, d] = numden(G); % n = numerator and d = denominator 

%roots:
eq = n == 0; % numerator must equal zero
roots_G = solve(eq, s);

% testing stability criteria:
if max(real(ev_G))>=0
    disp('Controlled system G is not asymptotically stable');
    stab_G = false;

elseif max(real(ev_G))<0
    disp('Controlled system G is asymptotically stable');
    stab_G = true;
end



end