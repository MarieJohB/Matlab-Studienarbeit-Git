% Helper function to interpret correlation strength
function strength = getCorrelationStrength(corr)
    if isnan(corr)
        strength = 'N/A';
        return;
    end
    
    absCorr = abs(corr);
    
    if absCorr >= 0.8
        strength = 'Very Strong';
    elseif absCorr >= 0.6
        strength = 'Strong';
    elseif absCorr >= 0.4
        strength = 'Moderate';
    elseif absCorr >= 0.2
        strength = 'Weak';
    else
        strength = 'Very Weak';
    end
    
    % Add direction
    if corr > 0
        strength = ['Positive, ' strength];
    elseif corr < 0
        strength = ['Negative, ' strength];
    end
end