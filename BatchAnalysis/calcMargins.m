function [Gm, Pm, Wcg, Wcp] = calcMargins(mag, phase, omega)
% Calculate gain and phase margins from Bode data

% Initialize outputs
Gm = NaN;
Pm = NaN;
Wcg = NaN;
Wcp = NaN;

% Convert phase to range -360 to 0
phase = mod(phase, -360);

% Find phase = -180° crossing (gain margin)
phase180idx = find(diff(sign(phase + 180)) ~= 0);

if ~isempty(phase180idx)
    % Get closest index
    phase180idx = phase180idx(1);
    
    % Linear interpolation to find exact frequency
    if phase180idx < length(phase)
        f1 = omega(phase180idx);
        f2 = omega(phase180idx + 1);
        p1 = phase(phase180idx) + 180;
        p2 = phase(phase180idx + 1) + 180;
        
        % Calculate exact frequency where phase = -180°
        Wcp = f1 + (f2 - f1) * (-p1) / (p2 - p1);
        
        % Find gain at this frequency
        g1 = mag(phase180idx);
        g2 = mag(phase180idx + 1);
        g_at_180 = g1 + (g2 - g1) * (-p1) / (p2 - p1);
        
        % Gain margin
        Gm = 1 / g_at_180;
    end
end

% Find magnitude = 1 (0 dB) crossing (phase margin)
mag0dBidx = find(diff(sign(20*log10(mag))) ~= 0);

if ~isempty(mag0dBidx)
    % Get closest index
    mag0dBidx = mag0dBidx(1);
    
    % Linear interpolation to find exact frequency
    if mag0dBidx < length(mag)
        f1 = omega(mag0dBidx);
        f2 = omega(mag0dBidx + 1);
        m1 = 20*log10(mag(mag0dBidx));
        m2 = 20*log10(mag(mag0dBidx + 1));
        
        % Calculate exact frequency where magnitude = 0 dB
        Wcg = f1 + (f2 - f1) * (-m1) / (m2 - m1);
        
        % Find phase at this frequency
        p1 = phase(mag0dBidx);
        p2 = phase(mag0dBidx + 1);
        p_at_0dB = p1 + (p2 - p1) * (-m1) / (m2 - m1);
        
        % Phase margin
        Pm = 180 + p_at_0dB;
    end
end
end