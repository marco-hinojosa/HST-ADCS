function [xhist, Phist] = mekf(x0, P0, W, V, rN, whist, yhist, dt)

xhist = zeros(7,size(yhist,2));
xhist(:,1) = x0;

Phist = zeros(6,6,size(yhist,2));
Phist(:,:,1) = P0;

for k = 1:(size(yhist,2)-1) 
    [xp, A] = prediction(xhist(:,k),whist(:,k),dt);
    P_p = A*Phist(:,:,k)*A' + W;   
    [yp, C] = measurement(xp(1:4),rN);   
    
    % Innovation
    z_q = quat2phi(qmult(qconj(xp(1:4)))*yhist(1:4,k+1)); % Computes ST Quaternion Innovation
    z_r = yhist(5:end,k+1) - yp(5:end); % Computes Body Vector Innovation
    z = [z_q; z_r];
    S = C*P_p*C' + V; 

    % Kalman Gain
    L = P_p*C'*inv(S); % Avoid Singularity
    
    % Update
    dx = L*z; % [dphi; dbeta], 6x1
    phi = dx(1:3);
    %dq = [1/2*phi ; 1 - 1/8*phi'*phi];
    dq = phi2quat(phi);
    dq = dq/norm(dq);
    
    % Quaternion and Bias Update
    xhist(1:4,k+1) = qmult(xp(1:4))*dq;
    xhist(5:7,k+1) = xp(5:7) + dx(4:6);
    
    % Covariance Update
    Phist(:,:,k+1) = (eye(6) - L*C)*P_p*(eye(6) - L*C)' + L*V*L';
end
end

