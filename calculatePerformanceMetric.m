function metricValue = calculatePerformanceMetric(sys, metric)
    % CALCULATEPERFORMANCEMETRIC - Calculate the specified performance metric
    % for a given system with improved error handling and gain margin calculation
    
    switch metric
        case 'Settling Time'
            try
                % Try to get step info
                info = stepinfo(sys);
                metricValue = info.SettlingTime;
            catch
                % If stepinfo fails, try a more robust approach
                metricValue = estimateSettlingTime(sys);
            end
            
        case 'Overshoot'
            try
                info = stepinfo(sys);
                metricValue = info.Overshoot;
            catch
                metricValue = NaN;
            end
            
        case 'Rise Time'
            try
                info = stepinfo(sys);
                metricValue = info.RiseTime;
            catch
                metricValue = NaN;
            end
            
        case 'Bandwidth'
            try
                metricValue = bandwidth(sys);
            catch
                % If bandwidth fails, use frequency response approach
                try
                    [mag, ~, w] = bode(sys, {0.01, 100});
                    mag = squeeze(mag);
                    idx = find(mag <= 0.707, 1, 'first');
                    if ~isempty(idx)
                        metricValue = w(idx);
                    else
                        metricValue = NaN;
                    end
                catch
                    metricValue = NaN;
                end
            end
            
        case 'Phase Margin'
            try
                [~, pm] = margin(sys);
                metricValue = pm;
            catch
                metricValue = NaN;
            end
            
        case 'Gain Margin'
            try
                % Improved gain margin calculation with better error handling
                [gm, ~] = margin(sys);
                
                % Check for infinite gain margin (stable system)
                if isinf(gm)
                    metricValue = 100; % Use a large value for stable systems
                else
                    % Convert to dB for consistency
                    metricValue = 20*log10(gm);
                    
                    % Handle unreasonably large values
                    if metricValue > 200
                        metricValue = 200;
                    end
                    
                    % Handle negative gain margins (indicates instability)
                    if metricValue < 0
                        metricValue = 0;
                    end
                end
            catch ME
                disp(['Gain Margin calculation error: ', ME.message]);
                metricValue = NaN;
            end
            
        otherwise
            metricValue = NaN;
    end
end