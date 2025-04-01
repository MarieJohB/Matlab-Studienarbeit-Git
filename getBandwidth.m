% Support function for bandwidth estimation
function w_bandwidth = getBandwidth(G, plantInfo)
    % GETBANDWIDTH Estimate bandwidth of plant
    
    try
        % Try standard bandwidth calculation for stable systems
        if ~plantInfo.isUnstable
            w_bandwidth = bandwidth(G);
            if ~isnan(w_bandwidth)
                return;
            end
        end
    catch
        % Continue to alternative methods if bandwidth calculation fails
    end
    
    % Alternative methods for bandwidth estimation
    try
        % Method 1: Use frequency response
        w = logspace(-3, 3, 500);
        [mag, ~] = bode(G, w);
        mag = squeeze(mag);
        
        % Find 0.7 (-3dB) crossing point
        try
            dc_gain = dcgain(G);
        catch
            % If DC gain calculation fails, use low-frequency gain
            dc_gain = mag(1);
        end
        
        threshold = 0.7 * abs(dc_gain);
        cross_idx = find(mag < threshold, 1);
        
        if ~isempty(cross_idx) && cross_idx > 1
            w_bandwidth = w(cross_idx);
            return;
        end
    catch
        % Continue to next method
    end
    
    % Method 2: Use pole/zero information
    try
        p = plantInfo.poles;
        
        if ~isempty(p)
            % Use dominant (slowest) pole for stable systems 
            % or fastest unstable pole for unstable systems
            if plantInfo.isUnstable
                unstable_poles = p(real(p) > 0);
                if ~isempty(unstable_poles)
                    [~, idx] = max(real(unstable_poles));
                    w_bandwidth = 5 * abs(unstable_poles(idx));  % Higher bandwidth needed for unstable pole
                    return;
                end
            else
                % For stable systems, use dominant poles
                stable_poles = p(real(p) < 0);
                if ~isempty(stable_poles)
                    [~, idx] = min(abs(real(stable_poles)));
                    w_bandwidth = 2 * abs(real(stable_poles(idx)));  % Rule of thumb
                    return;
                end
            end
        end
    catch
        % Continue to next method
    end
    
    % Method 3: Use FOPDT parameters if available
    if ~isnan(plantInfo.FOPDT.T) && plantInfo.FOPDT.T > 0
        w_bandwidth = 2.5 / plantInfo.FOPDT.T;  % Rule of thumb
        return;
    end
    
    % Default fallback value
    w_bandwidth = 1.0;
end