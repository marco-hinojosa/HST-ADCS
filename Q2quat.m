function q = Q2quat(Q)
% Converts Rotation Matrix Q into Quaternion q
% Input - Q Rotation Matrix 
% Output - q Quaternion
phi = unhat(logm(Q));
theta = norm(phi);
r = phi/theta;
q = [r*sin(theta/2); cos(theta/2)];