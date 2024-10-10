function is_sternlyproper = sternlyproper_tf(G)
% Function to check if a transfer function is sternly proper

syms s; % set 's' as symbolic variable

% Initialize output variables

[Num,Den] = tfdata(G,'v');

G_sys = poly2sym(Num,s)/poly2sym(Den,s); 

is_sternlyproper = false;

if limit(G_sys, s, inf) == 0
% Transfer function is sternly proper
is_sternlyproper = true;
end