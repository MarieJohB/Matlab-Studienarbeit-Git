function limitations = analyzeRHPZeroLimitations(plantInfo)
    % ANALYZERHPZEROLIMITATIONS Analyze fundamental limitations imposed by RHP zeros
    % Provides quantitative guidance on achievable performance with RHP zeros
    %
    % Input:
    %   plantInfo - Structure with plant information from analyzePlant()
    %
    % Output:
    %   limitations - Structure containing quantified performance limitations
    
    % Initialize limitations structure
    limitations = struct();
    limitations.hasRHPZeros = false;
    limitations.maxBandwidth = Inf;
    limitations.minRiseTime = 0;
    limitations.minSettlingTime = 0;
    limitations.maxIntegralGain = Inf;
    limitations.minPeakTime = 0;
    limitations.waterbed = [];
    limitations.recommendedSettings = [];
    limitations.description = '';
    
    % Get plant poles and zeros
    p = plantInfo.poles;
    z = plantInfo.zeros;
    
    % Extract RHP zeros (if any)
    rhp_zeros = z(real(z) > 0);
    
    if isempty(rhp_zeros)
        limitations.hasRHPZeros = false;
        limitations.description = 'Plant has no RHP zeros. No special limitations.';
        return;
    end
    
    % We have RHP zeros - need to analyze fundamental limitations
    limitations.hasRHPZeros = true;
    
    % Sort RHP zeros by real part
    [sorted_real, idx] = sort(real(rhp_zeros));
    sorted_rhp_zeros = rhp_zeros(idx);
    
    % Calculate the minimum RHP zero (closest to imaginary axis)
    min_real_rhp_zero = min(real(rhp_zeros));
    min_idx = find(real(rhp_zeros) == min_real_rhp_zero, 1);
    critical_zero = rhp_zeros(min_idx);
    
    % Calculate bandwidth limitation
    % Theoretical maximum bandwidth is approximately half the real part
    % of the rightmost RHP zero, but we apply a safety factor
    if imag(critical_zero) == 0
        % Real RHP zero
        limitations.maxBandwidth = min_real_rhp_zero * 0.5;
        limitations.safetyBandwidth = min_real_rhp_zero * 0.35; % More conservative
    else
        % Complex RHP zero - more restrictive
        omega_n = abs(critical_zero);
        zeta = real(critical_zero) / omega_n;
        limitations.maxBandwidth = min_real_rhp_zero * 0.5;
        limitations.safetyBandwidth = min_real_rhp_zero * 0.3; % More conservative for complex zeros
    end
    
    % Calculate minimum achievable rise time
    % According to theory, rise time is limited by RHP zeros
    limitations.minRiseTime = 1.5 / min_real_rhp_zero;
    
    % Calculate minimum settling time (more conservative than rise time)
    limitations.minSettlingTime = 3 / min_real_rhp_zero;
    
    % Calculate maximum recommended integral gain
    limitations.maxIntegralGain = min_real_rhp_zero * 0.2;
    
    % Minimum peak time (related to RHP zero locations)
    limitations.minPeakTime = 1 / min_real_rhp_zero;
    
    % Analyze waterbed effect - tradeoffs between sensitivity functions
    limitations.waterbed = analyzeSensitivityTradeoffs(critical_zero);
    
    % Get recommended controller settings for plants with RHP zeros
    limitations.recommendedSettings = getRecommendedSettings(sorted_rhp_zeros, p);
    
    % Create composite description
    desc = sprintf(['Plant has %d RHP zero(s) with minimum real part at %.3f rad/s.\n', ...
           'This imposes fundamental limitations on achievable performance:\n', ...
           '  • Maximum bandwidth: %.3f rad/s (conservative: %.3f rad/s)\n', ...
           '  • Minimum rise time: %.3f s\n', ...
           '  • Minimum settling time: %.3f s\n', ...
           '  • Maximum recommended integral gain: %.3f\n'], ...
           length(rhp_zeros), min_real_rhp_zero, ...
           limitations.maxBandwidth, limitations.safetyBandwidth, ...
           limitations.minRiseTime, limitations.minSettlingTime, ...
           limitations.maxIntegralGain);
    
    % Check for unstable poles to identify potential issues
    unstable_poles = p(real(p) > 0);
    
    if ~isempty(unstable_poles)
        % Sort unstable poles by real part
        [~, idx] = sort(real(unstable_poles));
        closest_unstable = unstable_poles(idx(1));
        
        % Check for potential design issues with RHP zeros and unstable poles
        if real(closest_unstable) > min_real_rhp_zero
            desc = [desc, sprintf('\nWARNING: Plant has unstable pole(s) with real part larger than RHP zero.\n')];
            desc = [desc, 'This combination creates severe limitations that may make control very difficult.\n'];
            desc = [desc, 'Consider redesigning the plant if possible.\n'];
        else
            % Calculate minimum difference between RHP zeros and unstable poles
            min_diff = Inf;
            for i = 1:length(rhp_zeros)
                for j = 1:length(unstable_poles)
                    diff = abs(rhp_zeros(i) - unstable_poles(j));
                    min_diff = min(min_diff, diff);
                end
            end
            
            if min_diff < 2
                desc = [desc, sprintf('\nWARNING: Plant has RHP zero(s) close to unstable pole(s) (distance: %.3f).\n', min_diff)];
                desc = [desc, 'This creates challenging control problems and may limit performance.\n'];
            end
        end
    end
    
    % Add information about multiple RHP zeros if present
    if length(rhp_zeros) > 1
        desc = [desc, sprintf('\nPlant has multiple RHP zeros, which compounds limitations.\n')];
        desc = [desc, 'Use a very conservative design approach.\n'];
    end
    
    limitations.description = desc;
end

function waterbed = analyzeSensitivityTradeoffs(critical_zero)
    % Analyze sensitivity function tradeoffs due to RHP zeros
    
    waterbed = struct();
    
    % Bode's sensitivity integral constraint due to RHP zero
    waterbed.integralConstraint = 'Bode sensitivity integral imposes that S(s) must have area under |S(jω)| > 0';
    
    % Get critical frequency of RHP zero
    z_real = real(critical_zero);
    
    % Calculate frequencies where sensitivity must be high
    if imag(critical_zero) ~= 0
        % Complex RHP zero
        waterbed.criticalFrequency = abs(critical_zero);
        waterbed.augmentedFrequency = waterbed.criticalFrequency * 2;
    else
        % Real RHP zero
        waterbed.criticalFrequency = z_real;
        waterbed.augmentedFrequency = z_real * 2;
    end
    
    % Calculate Poisson sensitivity integral constraint
    waterbed.poissonConstraint = pi * z_real;
    
    % Describe trade-offs
    waterbed.tradeoff = 'Reducing sensitivity at low frequencies will cause increased sensitivity at high frequencies';
    waterbed.recommendation = sprintf('Design for sensitivity peak < 6dB and accept poor high-frequency disturbance rejection above %.3f rad/s', waterbed.augmentedFrequency);
    
    return;
end

function settings = getRecommendedSettings(rhp_zeros, poles)
    % Get recommended controller settings for plants with RHP zeros
    
    settings = struct();
    
    % Find critical RHP zero (closest to imaginary axis)
    critical_zero = rhp_zeros(1);
    
    % Extract unstable poles (if any)
    unstable_poles = poles(real(poles) > 0);
    
    % Recommended bandwidth
    settings.bandwidth = real(critical_zero) * 0.3;
    
    % Recommended phase margin
    settings.phaseMargin = 60; % More conservative
    
    % Recommended damping ratio
    settings.damping = 0.8;   % More conservative
    
    % Recommended controller structure
    if isempty(unstable_poles)
        % For stable plants with RHP zeros
        settings.structure = 'PID with reduced integral action';
        settings.Kp = 'Moderate';
        settings.Ki = sprintf('Low (< %.3f)', real(critical_zero) * 0.2);
        settings.Kd = 'Moderate with heavy filtering';
        settings.specialNotes = 'Use 2nd order filtering on derivative term';
    else
        % For unstable plants with RHP zeros
        settings.structure = 'Compensator with reduced integral action';
        settings.Kp = 'Low';
        settings.Ki = sprintf('Very low (< %.3f)', real(critical_zero) * 0.1);
        settings.Kd = 'Higher with heavy filtering';
        settings.specialNotes = 'Consider plant redesign if possible, as this combination is fundamentally difficult';
    end
    
    % Recommended robustness level
    settings.robustness = 'High';
    
    % Recommendations for different control goals
    settings.trackingRecommendation = 'Accept slow response times and prioritize overshoot reduction';
    settings.disturbanceRecommendation = 'Accept limited low-frequency rejection and avoid aggressive integral action';
    settings.robustnessRecommendation = 'Use conservative gain margins (>10dB) and phase margins (>60 deg)';
    
    return;
end