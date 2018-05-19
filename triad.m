function [q,Q] = triad(rN,rB)
% TRIAD Algorithm for Static State Estimation
% Inputs:
%      rN1 - First Inertial-Frame Reference Position Vector
%      rN2 - Second Inertial-Frame Reference Position Vector
%      rB1 - First Body-Frame Measured Position Vector
%      rB2 - Second Body-Frame Measured Position Vector
%
% Outputs:
%   Q - 3x3 Rotation Matrix
%   q - 4x1 Attitude Quaternion

rN1 = rN(:,1); rN2 = rN(:,2); 
rB1 = rB(:,1); rB2 = rB(:,2); 
tn1 = rN1;
tn2 = cross(rN1,rN2)/norm(cross(rN1,rN2));

Tn = [tn1 tn2 cross(tn1,tn2)/norm(cross(tn1,tn2))];

tb1 = rB1;
tb2 = cross(rB1,rB2)/norm(cross(rB1,rB2));

Tb = [tb1 tb2 cross(tb1,tb2)/norm(cross(tb1,tb2))];

% Compute Q and q
Q = Tn*Tb';
q = Q2quat(Q);
end