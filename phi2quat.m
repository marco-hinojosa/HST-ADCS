function q = phi2quat(phi)
% Converts Axis-Angle Vector phi into Quaternion q
% Input - phi Axis Angle 
% Output - q Quaternion Matrix
if norm(phi) > 10*pi/180 
    theta = norm(phi);
    r = phi/theta;

    q(1:3) = r.*sin(theta/2);
    q(4) = cos(theta/2);
    q = q';
else % Small Angle Approximation
    q = [phi./2 ; 1 - 1/8*(phi'*phi)];
    q = q/norm(q);
end
end