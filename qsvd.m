function [q_svd,Q_svd] = qsvd(rN,rB,stddev)
% SVD Algorithm for Static State Estimation
% Inputs:
%      rN1 - First Inertial-Frame Reference Position Vector
%      rN2 - Second Inertial-Frame Reference Position Vector
%      rB1 - First Body-Frame Measured Position Vector
%      rB2 - Second Body-Frame Measured Position Vector
%
% Outputs:
%   Q - 3x3 Rotation Matrix
%   q - 4x1 Attitude Quaternion

% Generate Attitude Profile Matrix
B = zeros(3);

for ii=1:size(rB,2)
    B = B + rB(:,ii)*rN(:,ii)'/stddev;
end

[U,S,V] = svd(B');
Q_svd = U*diag([1,1,det(U)*det(V)])*V';
q_svd = Q2quat(Q_svd);

% Generate Rotation Error
% Q_error = Q_true'*Q_NB_SVD;
% e = unhat(logm(Q_error));
% th_e = 180/pi*norm(e);