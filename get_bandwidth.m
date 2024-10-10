function bandwidth_T = get_bandwidth(T)
%Calculate the Bandwidth of T

    try
        bandwidth_T = bandwidth(T);
        if isnan(bandwidth_T)
            disp('warning: check if model has infinite DC gain. Searching for -3 dB bandwidth instead.');
            bandwidth_T = NaN; % Manuelles Setzen von NaN als Indikator
        end
        
    catch 
        disp('error trying to find bandwidth');
        bandwidth_T = NaN; 
    end

    % if no bandwidth could be found, check for -3dB point
    if isnan(bandwidth_T)
        [mag, ~, freq] = bode(T);
        mag = squeeze(mag);
        bandwidth_idx = find(20*log10(mag) <= -3, 1, 'first');
        if isempty(bandwidth_idx)
            disp('No frequency crosses -3 dB. Bandwidth cannot be defined.');
            bandwidth_T = NaN;
        else
            bandwidth_T = freq(bandwidth_idx);
        end
    end


end