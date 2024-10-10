function get_bode(L)
% Function to plot the Bode diagram of a transfer function and calculate gain margin, phase margin, and crossover frequency

% Check if the input is a valid transfer function
if ~isa(L, 'tf')
  error('Input must be a transfer function object');
end
  
% Calculate margins
% Gm = gain margin 
% Pm = phase margin
% wcg = gain crossover frequency
% wcp = phase crossover frequency
[Gm, Pm, Wcg, Wcp] = margin(L);
    
margin(L);


% Output the calculated values
fprintf('Gain Margin: %.2f dB\n', 20*log10(Gm));
fprintf('Phase Margin: %.2f°\n', Pm);
fprintf('Crossover Frequency: %.2f rad/s\n', Wcg);

end