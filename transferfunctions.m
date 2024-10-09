function [T, S, L] = transferfunctions(G, K)

% function to calculate different transferfunctions 
% inputs: controlled system G(s) and controller K(s)
% outputs: complementary sensitivity T, sensitivity S, open loop L

T = feedback(G*K,1);
S = feedback(1,G*K);
L = G*K;
end