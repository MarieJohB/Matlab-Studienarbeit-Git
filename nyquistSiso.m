function [lReal, lImag, omega] = nyquistSiso(L, omega)
    % [ lReal, lImag, omega ] = nyquistSiso(L, omega)
    % Computes and displays the Nyquist diagram of a SISO system L 
    % for positive angular frequencies omega

    if nargin == 1
        [lReal, lImag, omega] = nyquist(L); 
    elseif nargin == 2
        [lReal, lImag] = nyquist(L, omega); 
    end
    lReal = squeeze(lReal); 
    lImag = squeeze(lImag); 

end

