function unit = getMetricUnit(metric)
    switch metric
        case 'Settling Time'
            unit = ' s';
        case 'Overshoot'
            unit = ' %';
        case 'Rise Time'
            unit = ' s';
        case 'Bandwidth'
            unit = ' rad/s';
        case 'Phase Margin'
            unit = ' degrees';
        case 'Gain Margin'
            unit = ' dB';
        otherwise
            unit = '';
    end
end