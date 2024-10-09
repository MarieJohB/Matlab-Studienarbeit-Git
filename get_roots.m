function [results_roots] = get_roots(A)
% input:  tranfer function
% output: roots_G: roots of the function 
syms s; % set 's' as sybolic variable

% calculation of roots (Nullstellen) of function
 
[n, d] = tfdata(A); % n = numerator and d = denominator 

%roots:
num = cell2mat(n);

results_roots = roots(num);

end