function [bw, cf] = get_bode(A)
% input: transfer fuction or element to by analysed
% output: 
% bode plot
% bandwidth 
% crossover frequency 

bw = bandwidth(A);

% bode plot
figure;
bode(A);
title('Bode Diagramm');


end