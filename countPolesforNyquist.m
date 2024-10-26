function [m_0, n_0, deltaphi] = countPolesforNyquist(L)
    % Get poles of the transfer function
    poles = pole(L);
    
    % Initialize counters
    m_0 = 0; % Poles with real part > 0
    n_0 = 0; % Poles on the imaginary axis
    
    % Iterate through poles to count
    for k = 1:length(poles)
        if real(poles(k)) > 0
            m_0 = m_0 + 1;
        end
        if real(poles(k)) == 0 && imag(poles(k)) ~= 0 | real(poles(k)) == 0 && imag(poles(k)) == 0
            n_0 = n_0 + 1;
        end
    end
    
    % Calculate deltaphi
    deltaphi = m_0 * 180 + n_0 * 90;
end
