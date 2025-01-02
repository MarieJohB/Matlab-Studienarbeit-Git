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

    % graphical output
    plot(lReal, lImag, 'b', 'LineWidth', 2);
    grid on;
    hold on;
    plot(lReal(1), lImag(1), 'oc', 'LineWidth', 2);
    plot(lReal(end), lImag(end), 'xc', 'LineWidth', 2);
    plot(-1, 0, '+r', 'LineWidth', 2);
    hold off;
    xlabel('Real Axis');
    ylabel('Imaginary Axis');
    title('Nyquist Diagram for \omega=0...+\infty, o=L(j0), x=L(j\infty), +=critical point');
end
