function is_proper = proper_tf(G)

    % Get numerator and denominator order
    [num, den] = tfdata(G, 'v');
    numeratorOrder = length(num) - 1;
    denominatorOrder = length(den) - 1;
    
    % Check if the transfer function is proper
    is_proper = numeratorOrder <= denominatorOrder;
end
