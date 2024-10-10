function get_bode(G)
% Function to plot the Bode diagram of a transfer function and calculate gain margin, phase margin, and crossover frequency

% Check if the input is a valid transfer function
if ~isa(G, 'tf')
  error('Input must be a transfer function object');
end
  
% Calculate margins
[Gm, Pm, Wcg, Wcp] = margin(G);
    
margin(G);

% bode plot
figure;
bode(A);
title('Bode Diagramm');


% Output the calculated values
fprintf('Gain Margin: %.2f dB\n', 20*log10(Gm));
fprintf('Phase Margin: %.2f°\n', Pm);
fprintf('Crossover Frequency: %.2f rad/s\n', Wcg);

end