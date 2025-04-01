function rhpLimits = analyzeRHPZeroLimitations(plantInfo)
    % ANALYZERHPZEROLIMITATIONS Calculate performance limitations
    % due to RHP zeros
    
    rhpLimits = struct();
    
    % Extract RHP zeros
    rhpZeros = plantInfo.zeros(real(plantInfo.zeros) > 0);
    
    if isempty(rhpZeros)
        return;
    end
    
    % Find the minimum real RHP zero
    realRHPZeros = rhpZeros(imag(rhpZeros) == 0);
    
    if ~isempty(realRHPZeros)
        minRealRHPZero = min(realRHPZeros);
        
        % Calculate bandwidth limitation
        % Rule of thumb: bandwidth < 0.5 * min RHP zero
        rhpLimits.maxBandwidth = 0.5 * minRealRHPZero;
        
        % Calculate rise time limitation
        % Rule of thumb: rise time > 1.5 / min RHP zero
        rhpLimits.minRiseTime = 1.5 / minRealRHPZero;
        
        % Calculate sensitivity peak bound
        % Approximate Ms > 1 + 2/|z|, where |z| is the distance to the origin
        rhpLimits.minSensitivityPeak = 1 + 2/minRealRHPZero;
    else
        % For complex RHP zeros, just use the real part
        minRealPart = min(real(rhpZeros));
        
        rhpLimits.maxBandwidth = 0.5 * minRealPart;
        rhpLimits.minRiseTime = 1.5 / minRealPart;
        rhpLimits.minSensitivityPeak = 1 + 2/minRealPart;
    end
    
    % Check for PIP (Parity Interlacing Property) violations
    if plantInfo.isUnstable
        rhpPoles = plantInfo.poles(real(plantInfo.poles) > 0);
        
        % Sort RHP zeros and poles
        sortedRHPZeros = sort(real(rhpZeros));
        sortedRHPPoles = sort(real(rhpPoles));
        
        % Check interlacing property
        rhpLimits.pipViolation = false;
        
        if length(sortedRHPZeros) >= 1 && length(sortedRHPPoles) >= 1
            % Check if each zero has at least one pole to its right
            for i = 1:length(sortedRHPZeros)
                polesRight = sum(sortedRHPPoles > sortedRHPZeros(i));
                if polesRight < i
                    rhpLimits.pipViolation = true;
                    break;
                end
            end
        end
        
        if rhpLimits.pipViolation
            rhpLimits.pipViolationNote = 'PIP violation detected. This system may be difficult to control with standard methods. Consider using Youla-Kucera or Pre-stabilization.';
        else
            rhpLimits.pipViolationNote = '';
        end
    else
        rhpLimits.pipViolation = false;
        rhpLimits.pipViolationNote = '';
    end
end