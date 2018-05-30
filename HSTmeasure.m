function [y,beta] = HSTmeasure(x,V,flag,rN)
% HSTmeasure: Creates sensor measurements by adding noise to true states.
%
% Inputs:
%      x - Ttime History of True States
%      V - Sensor Noise Covariance
%      flag - Sensor Type
%
% Outputs:
%   y - Noisy Sensor Measurements
%   beta - Bias History (for Gyros)

if isempty(rN)
else
    rN1 = rN(:,1);
    rN2 = rN(:,2);
end
if flag == 1 % Magnetometers Body Vectors
    for ii = 1:length(x)
        q_stat = x(ii,4:7)';
        Q = quat2Q(q_stat);
        
        rB1(:,ii) = (Q'*rN1) + chol(V)*randn(3,1);
        rB1(:,ii) = rB1(:,ii)/norm(rB1(:,ii));
        rB2(:,ii) = (Q'*rN2) + chol(V)*randn(3,1);
        rB2(:,ii) = rB2(:,ii)/norm(rB2(:,ii));
        
        y(1:3,ii) = rB1(:,ii);
        y(4:6,ii) = rB2(:,ii);
    end
    beta = [];
elseif flag == 2 % Star Tracker (Quaternion)
    for ii = 1:length(x)
        phi = chol(V)*randn(3,1); % Axis-Angle Noise
        dq = phi2quat(phi); % Convert to Quaternion Noise
        dq = dq/norm(dq); % Re-Normalize
        y(:,ii) = qmult(x(ii,4:7)')*dq; % "Add" Noise by Multiplying to State Quaternion
    end
    beta = [];
elseif flag == 3 % Gyroscope (Vector)
    beta(:,1) = [0;0;0];
    W_rnd = V(1:3,1:3); % Rate Noise Density Covariance
    W_arw = V(4:6,4:6); % Angle Random Walk Covariance
    for ii = 1:length(x)
        om_true = x(ii,1:3)';
        beta(:,ii+1) = beta(:,ii) + chol(W_arw)*randn(3,1);
        om_noise(:,ii) = om_true + beta(:,ii) + chol(W_rnd)*randn(3,1);
    end
    y = om_noise;
    beta = beta(:,1:end-1);
end
end