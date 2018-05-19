function [q_dav,Q_dav] = qdaven(rN,rB,stddev)
% Davenport q-Algorithm for Static State Estimation
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
Z = zeros(3,1);

for ii=1:size(rB,2)
    B = B + rB(:,ii)*rN(:,ii)'/stddev;
    Z = Z + cross(rB(:,ii),rN(:,ii))/stddev;
end

K = [B + B' - trace(B)*eye(3) Z;Z' trace(B)];
[eigvec, eigval] = eig(K);
[~, idx] = max(diag(eigval));

q_dav = eigvec(:,idx)/norm(eigvec(:,idx));
Q_dav = quat2Q(q_dav);

% Generate Rotation Error
%Q_error = Q_true'*Q_dav;
%e = unhat(logm(Q_error));
%th_e = 180/pi*norm(e)