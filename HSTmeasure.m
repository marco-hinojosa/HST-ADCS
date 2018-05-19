function [y,v] = HSTmeasure(x,V,flag,rN)
% HSTmeasure: Creates sensor measurements by adding noise to true states.
%
% Inputs:
%      x - time history of true states
%      V - sensor noise covariance
%      flag - sensor type
%
% Outputs:
%   y - Noisy sensor measurements
%   v - Time history of sensor noise

if isempty(rN)
else
    rN1 = rN(:,1);
    rN2 = rN(:,2);
end
if flag == 1 % Magnetometers Body Vectors
    for ii = 1:length(x);
        q_stat = x(ii,4:7)';
        Q = quat2Q(q_stat);
        v(:,ii) = chol(V)*randn(3,1);
        
        rB1(:,ii) = (Q'*rN1) + v(:,ii);
        rB1(:,ii) = rB1(:,ii)/norm(rB1(:,ii));
        rB2(:,ii) = (Q'*rN2) + v(:,ii);
        rB2(:,ii) = rB2(:,ii)/norm(rB2(:,ii));
        
        y(:,ii,1) = rB1(:,ii);
        y(:,ii,2) = rB2(:,ii);
    end
elseif flag == 2; % Star Tracker (Quaternion)
    for ii = 1:length(x);
        phi = chol(V)*randn(3,1); % Axis-Angle Noise
        q_noise = [phi/2; 1-(1/8)*phi'*phi]; % Convert to Quaternion Noise
        q_noise = q_noise/norm(q_noise); % Renormalize
        v(:,ii) = q_noise; % Store Time History
        y(:,ii) = qmult(x(ii,4:7)')*v(:,ii); % "Add" Noise by Multiplying to State Quaternion
    end
elseif flag == 3; % Gyroscope (Vector)
    beta(:,1) = [0;0;0];
    V_rnd = V(:,:,1); % Rate Noise Density Covariance
    V_arw = V(:,:,2); % Angle Random Walk Covariance
    for ii = 1:length(x);
        om_true = x(ii,1:3)';
        v(:,ii) = chol(V_rnd)*randn(3,1);
        beta(:,ii+1) = beta(:,ii) + chol(V_arw)*randn(3,1);
        om_noise(:,ii) = om_true + beta(:,ii) + v(:,ii);
    end
    y = [om_noise;beta(:,1:end-1)];
end
end