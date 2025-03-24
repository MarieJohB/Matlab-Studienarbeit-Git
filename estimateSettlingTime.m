function settlingTime = estimateSettlingTime(sys)
    % ESTIMATESETTLINGTIME - More robust estimation of settling time
    % when standard stepinfo approach fails
    
    try
        % Simulate step response
        t = linspace(0, 100, 1000);
        [y, t] = step(sys, t);
        
        if isempty(y) || all(isnan(y))
            settlingTime = NaN;
            return;
        end
        
        % Get final value (using last 10% of response)
        finalIdx = round(0.9 * length(y)):length(y);
        if isempty(finalIdx)
            finalValue = y(end);
        else
            finalValue = mean(y(finalIdx));
        end
        
        % Find when response stays within 2% of final value
        tolerance = 0.02 * abs(finalValue);
        withinBand = abs(y - finalValue) <= tolerance;
        
        % Find the last time point where the response leaves the band
        if any(withinBand)
            lastOutsideBandIdx = find(~withinBand, 1, 'last');
            if isempty(lastOutsideBandIdx)
                settlingTime = t(1); % Always within band
            else
                settlingTime = t(lastOutsideBandIdx);
            end
        else
            settlingTime = NaN; % Never settles
        end
    catch
        settlingTime = NaN;
    end
end