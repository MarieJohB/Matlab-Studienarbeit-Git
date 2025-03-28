function cmap = redblue(m)
    % REDBLUE Creates a red-blue colormap centered at zero
    % This function creates a colormap with red for negative values,
    % white for zero, and blue for positive values.
    %
    % Parameters:
    %   m - Number of color entries (default: 64)
    %
    % Returns:
    %   cmap - m-by-3 colormap matrix
    
    if nargin < 1
       m = 64;
    end
    
    % Calculate midpoint
    mid = floor(m/2);
    
    % Create red part (negative values)
    r = ones(mid, 1);
    g = linspace(0, 1, mid)';
    b = linspace(0, 1, mid)';
    
    % Create blue part (positive values) 
    r2 = linspace(1, 0, m-mid)';
    g2 = linspace(1, 0, m-mid)';
    b2 = ones(m-mid, 1);
    
    % Combine to create final colormap
    cmap = [r g b; r2 g2 b2];
end