function [ev_G, stab_G] = stability_eigenvalues(G)
% caluculation of eigenvalues
% checking for stability of the controlled systems based on calculation of
% eigenvalues / poles  N_G = 0

% input: G (controlled system)
% output: 
% ev_G: eigenvalues of G
% stab_G: status whether G is stable or not


syms s; % set 's' as sybolic variable

% calculation of eigenvalues / poles of G
ev_G = eig(G);

% testing stability criteria:
if max(real(ev_G))>=0
    disp('System is not asymptotically stable');
    stab_G = false;

elseif max(real(ev_G))<0
    disp('System is asymptotically stable');
    stab_G = true;
end

end