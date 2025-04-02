% Generate unique desired poles based on pole type and parameters
function poles = generateUniqueDesiredPoles(poleType, n_states, params)
    % Initialize poles array
    poles = zeros(1, n_states);
    
    switch poleType
        case 'Real Poles'
            % Generate evenly spaced real poles with bandwidth as the base
            bandwidth = params.bandwidth;
            
            % Get spacing parameter or use default
            if isfield(params, 'spacing')
                spacing = params.spacing;
            else
                spacing = 0.5; % Default spacing
            end
            
            % Generate real poles with spacing
            for i = 1:n_states
                poles(i) = -bandwidth - (i-1) * spacing;
            end
            
        case 'Complex Poles (2nd order)'
            % Generate complex conjugate pairs based on bandwidth and damping
            bandwidth = params.bandwidth;
            
            if isfield(params, 'damping')
                damping = params.damping;
            else
                damping = 0.7; % Default damping
            end
            
            % Calculate natural frequency from bandwidth
            wn = bandwidth / (damping + sqrt(damping^2 + 1)); % Approximation
            
            % Number of pairs needed
            n_pairs = floor(n_states / 2);
            
            % Generate complex conjugate pairs
            for i = 1:n_pairs
                % Calculate current pair with slight variation in natural frequency
                % to ensure uniqueness for SISO systems
                current_wn = wn * (1 + 0.1 * (i-1));
                
                % Calculate complex pole pair
                real_part = -damping * current_wn;
                imag_part = current_wn * sqrt(1 - damping^2);
                
                % Add the pair
                poles(2*i-1) = complex(real_part, imag_part);
                poles(2*i) = complex(real_part, -imag_part);
            end
            
            % If odd number of states, add one real pole
            if mod(n_states, 2) == 1
                poles(n_states) = -bandwidth;
            end
            
        case 'Combined Poles'
            % A mix of real and complex poles
            bandwidth = params.bandwidth;
            
            if isfield(params, 'damping')
                damping = params.damping;
            else
                damping = 0.7; % Default damping
            end
            
            % Calculate natural frequency
            wn = bandwidth;
            
            % Roughly half real, half complex
            n_complex_pairs = floor(n_states / 4);
            n_real = n_states - 2 * n_complex_pairs;
            
            % Generate complex pairs
            for i = 1:n_complex_pairs
                current_wn = wn * (1 + 0.1 * (i-1));
                real_part = -damping * current_wn;
                imag_part = current_wn * sqrt(1 - damping^2);
                
                poles(2*i-1) = complex(real_part, imag_part);
                poles(2*i) = complex(real_part, -imag_part);
            end
            
            % Generate real poles
            for i = 1:n_real
                poles(2*n_complex_pairs + i) = -bandwidth * (1 + 0.2 * (i-1));
            end
    end
    
    % Ensure poles are unique for SISO systems
    if length(unique(poles)) < length(poles)
        % Add small variations to make poles unique
        for i = 1:n_states
            poles(i) = poles(i) * (1 + 0.01 * i);
        end
    end
end