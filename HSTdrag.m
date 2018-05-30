function [Drag,Torque] = HSTdrag(rvec,vvec,nQb,S,n,c)
rvec = rvec*1000; % [m]
vvec = vvec*1000; % [m/s]
rho = 2.72e-12; % [kg/m^3]
omE = (7.2921e-05)*[0; 0; 1]; % [rad/s]
Vn = (vvec + cross(omE,rvec)); % [m/s]
Vb = nQb'*Vn; % [m/s]
Cd = 2.2;

% Compute Atmospheric Drag and Torque (HST-Specific Values)
D = zeros(size(c));
Drag = zeros(3,1);
T = zeros(size(c));
Torque = zeros(3,1);
for ii = 1:size(c,2)
    D(:,ii) = .5*rho*Cd*norm(Vb)*Vb*S(ii)*max(dot(Vb,n(:,ii))/norm(Vb),0);
    T(:,ii) = cross(c(:,ii),D(:,ii));
    
    Drag = Drag + D(:,ii);
    Torque = Torque + T(:,ii);
end
