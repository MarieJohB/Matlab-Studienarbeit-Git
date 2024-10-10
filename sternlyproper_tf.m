function is_sternlyproper = sternlyproper_tf(G)
% Function to check if a transfer function is sternly proper

syms s; % set 's' as symbolic variable

% Initialize output variables

is_sternlyproper = false;

if limit(G, s, inf) == 0
% Transfer function is sternly proper
is_sternlyproper = true;
end