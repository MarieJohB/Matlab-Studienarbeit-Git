function is_sternlyproper = sternlyproper_tf(G)
% Function to check if a transfer function is sternly proper

syms s; % set 's' as symbolic variable

% Initialize output variables

[num,den] = tfdata(G,'v');

% creating symbolic transferfunction 
G_sys = poly2sym(num,s)/poly2sym(den,s); 



if limit(G_sys, s, inf) == 0
% Transfer function is sternly proper
is_sternlyproper = true;
%disp('Function is sternly proper');
else
is_sternlyproper = false;
%disp('Function is NOT sternly proper');
end 

end