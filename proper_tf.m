function is_proper = proper_tf(G)
% Function to check if a transfer function is proper

syms s; % set 's' as symbolic variable

[Num,Den] = tfdata(G,'v');

G_sys = poly2sym(Num,s)/poly2sym(Den,s); 

if limit(G_sys, s, inf) < inf
    disp('Transferfunction is proper');
    is_proper = true;
else
    disp('Transferfunction is not proper');
    is_proper = false;
end