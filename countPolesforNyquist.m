function [m_0, n_0, deltaphi] = countPolesforNyquist(L)
    % COUNTPOLESFORNYQUIST
    % -------------------
    % Counts the number of poles in the RHP and on the imaginary axis
    % to determine the required delta phi for Nyquist stability assessment
    
    % Get poles of the transfer function with improved numerical threshold
    poles = pole(L);
    
    % Define a small threshold for numerical stability
    epsilon = 1e-10;
    
    % Initialize counters
    m_0 = 0; % Poles with real part > 0
    n_0 = 0; % Poles on the imaginary axis
    
    % Iterate through poles to count with improved numerical handling
    for k = 1:length(poles)
        % Check for RHP poles (positive real part)
        if real(poles(k)) > epsilon
            m_0 = m_0 + 1;
        end
        
        % Check for imaginary axis poles (real part ≈ 0)
        if abs(real(poles(k))) <= epsilon
            n_0 = n_0 + 1;
        end
    end
    
    % Calculate required deltaphi in degrees
    deltaphi = m_0 * 180 + n_0 * 90;
    
    % Display information for debugging
    disp(['Open-loop poles: ', num2str(poles')]);
    disp(['m_0 (RHP poles): ', num2str(m_0)]);
    disp(['n_0 (imaginary axis poles): ', num2str(n_0)]);
    disp(['Required Δφ: ', num2str(deltaphi), '°']);
end