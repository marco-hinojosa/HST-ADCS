function [x_dot] = ODEstep(t,x,u,A,B,q,params)
% ODEstep: A single integration step for linearized HST attitude.
%
% Inputs:
%       t - time [s]
%       x0 - initial state vector [omega; q; r; v; speed]
%       mu - central body gravitational parameters [km^3/s^2]
%       J - spacecraft inertia matrix
%       params - environmental disturbance parameters
%
% Outputs:
%       x_dot - linear equation containing ODEs for all states

% Orbit Info
mu = 3.986e5; % Earth Standard Gravitational Parameter (km^3/s^2)
rvec = x(7:9);
vvec = x(10:12);

% Attitude Info
J11 = 28525.53;
J22 = 174815.86;
J33 = 181630.81;
J = diag([J11 J22 J33]); % Principal Axes of Inertia
Q = quat2Q(q); % nQb

% Gravity Gradient
Tau_g = cross(3*mu/(dot(rvec,rvec)^(5/2))*rvec, Q*J/(1000^2)*rvec); % [N-km]
Tau_g = Q'*Tau_g*1000; % Rotate into Body Frame and Convert Units, [N-m]

% Aero Drag
S = params(1,:);
n = params(2:4,:);
c = params(5:7,:);
[D,Tau_D] = HSTdrag(rvec,vvec,Q,S,n,c);
D = Q*D/1000; % Rotate Drag into Inertial Frame, [kg-km/s^2]
M = 11100; % HST Mass, [kg]

% Differential Equations for Linearized Phi/Omega, and r and v
B_tau = [zeros(3); inv(J)];
x_dot = [A*x(1:6) + [B B_tau]*[u; (Tau_g + Tau_D)]; vvec; -(mu/(norm(rvec)^3))*rvec - (D/M)];
end