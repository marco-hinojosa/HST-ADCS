%% AA279C ADCS Report
% Marco Hinojosa
% 06181747
close all,clear all,clc

%% Model HST Orbital Dynamics
% Reference Orbital Elements and Parameters
a = 6917.5; % Orbit Semi-Major Axis (km)
e = 0.000287; % Eccentricity
i = 28.47; % Inclination (deg)
RAAN = 176.23; % Right Ascension of Ascending Node (deg)
w = 82.61; % Argument of Perigee (deg)
anom = 319.41; % Mean Anomaly (deg)
mu = 3.986e5; % Earth Standard Gravitational Parameter (km^3/s^2)
T = 2*pi*sqrt((a^3)/mu); % Orbital Period (s)

% Convert to Earth-Centered Inertial Frame Coordinates
[r_eci, v_eci] = OE2ECI(a, e, i, RAAN, w, anom, mu);
x0 = [r_eci;v_eci];

% ODE45 Solver
tol = 1e-6;
opts  = odeset('reltol', tol, 'abstol', tol);
tspan = [0:10:3*T];
[t, x] = ode45(@(t,x)HSTorbit(mu,x), tspan, x0, opts);

% Plot Orbit
figure, hold on
plot3(x(:,1),x(:,2),x(:,3),'r','LineWidth',2)
axis equal
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')
title('Unperturbed HST Orbit')
earthPlot(1)
hold off
