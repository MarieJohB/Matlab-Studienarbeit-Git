function [proper, sternlyproper] = proper_tf(G)
% Function to check if a transfer function is strictly proper or proper

syms s; % set 's' as symbolic variable

% Initialize output variables
proper = false;
sternlyproper = false;

    if limit(G, s, inf) == 0
        % Transfer function is strictly proper
        sternlyproper = true;
        proper = true;
    elseif limit(G, s, inf) < inf
        % Transfer function is proper
        sternlyproper = false;
        proper = true;
    end
end