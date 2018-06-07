function [K,S] = tvlqr(A,B,Q,R,S)
%TVLQR Time-Varying Linear-Quadratic Regulator design from continuous 
%      cost function.

K = inv(R + B'*S*B)*B'*S*A;
S = Q + K'*R*K + (A - B*K)'*S*(A - B*K);