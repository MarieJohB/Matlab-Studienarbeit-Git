function is_proper = proper_tf(G)
% Function to check if a transfer function is proper

syms s; % set 's' as symbolic variable

% Initialize output variables
is_proper = false;

if limit(G, s, inf) < inf
  is_proper = true;
end