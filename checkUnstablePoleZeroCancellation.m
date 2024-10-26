function isCancellation = checkUnstablePoleZeroCancellation(G, K)
   % Calculate poles and zeroes
    G_poles = pole(G);
    G_zeroes = zero(G);
    K_poles = pole(K);
    K_zeroes = zero(K);

    % Initialize cancellation flag and details
    isCancellation = false;
    cancelDetails = [];
    
    % Check for cancellations between G and K
    for i = 1:length(G_zeroes)
        if any(abs(G_zeroes(i) - K_poles) < 1e-6)
            isCancellation = true;
            cancelDetails = [cancelDetails; G_zeroes(i)];
        end
    end
    
    for i = 1:length(K_zeroes)
        if any(abs(K_zeroes(i) - G_poles) < 1e-6)
            isCancellation = true;
            cancelDetails = [cancelDetails; K_zeroes(i)];
        end
    end
end