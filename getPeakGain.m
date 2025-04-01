function [peakgain, wpeak, w] = getPeakGain(sys)
    % Berechnet den maximalen Amplitudengang eines Systems
    w = logspace(-3, 3, 1000);
    [mag, ~] = bode(sys, w);
    mag = squeeze(mag);
    [peakgain, idx] = max(mag);
    wpeak = w(idx);
end