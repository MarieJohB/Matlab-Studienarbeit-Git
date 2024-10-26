function is_sternlyproper = sternlyproper_tf(G)
    % Get numerator and denominator order
    [num, den] = tfdata(G, 'v');
    numeratorOrder = length(num) - 1;
    denominatorOrder = length(den) - 1;

    % Check if the transfer function is sternly proper
    is_sternlyproper = numeratorOrder < denominatorOrder;
end