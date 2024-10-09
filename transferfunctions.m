function [T, S, L, GS, KS] = transferfunctions(G, K)

% function to calculate different transferfunctions 
% inputs: controlled system G(s) and controller K(s)
% outputs: transferfunctions of closed lop: complementary sensitivity T, sensitivity S, GS, KS
% outputs: transferfunction of open loop L

T = feedback(G*K,1);
S = feedback(1,G*K);

GS = minreal(G*S); % evtl aufpassen wegen Wegkuerzung?
KS = minreal(K*S);

L = G*K;
end