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

%% Model HST Rigid Body Dynamics (All Calcs Done in Principal Frame)
% Principal Axes of Inertia
J11 = 28525.53;
J22 = 174815.86;
J33 = 181630.81;
Jvec = [J11 J22 J33];
J = diag(Jvec);

% Initialize Parameters and Conditions
om_max = 10; % RPM
om_vec1 = [om_max;0;0];
om_vec1(1) = sqrt((om_max*J33)^2 - (J22*om_vec1(2))^2 - (J33*om_vec1(3))^2)/J11;
om_vec2 = [0;om_max;0];
om_vec2(2) = sqrt((om_max*J33)^2 - (J11*om_vec2(1))^2 - (J33*om_vec2(3))^2)/J22;
om_vec3 = [0;0;om_max];

h = norm(J*om_vec3);

om01 = [10;1;1];
om01(1) = sqrt((om_max*J33)^2 - (J22*om01(2))^2 - (J33*om01(3))^2)/J11;
om02 = [1;10;1];
om02(2) = sqrt((om_max*J33)^2 - (J11*om02(1))^2 - (J33*om02(3))^2)/J22;
om03 = [5;1;10];
om03(3) = sqrt((om_max*J33)^2 - (J11*om03(1))^2 - (J22*om03(2))^2)/J33;

% Integrate Trajectories
tol = 1e-6;
opts  = odeset('reltol', tol, 'abstol', tol);
tspan = [0:0.1:20];
[t, om1] = ode45(@(t,om)HSTspin(t,om,J), tspan, om01, opts);
[t, om2] = ode45(@(t,om)HSTspin(t,om,J), tspan, om02, opts);
[t, om3] = ode45(@(t,om)HSTspin(t,om,J), tspan, om03, opts);

[X,Y,Z] = sphere(50);
CO(:,:,1) = 0.6.*ones(size(X,1));
CO(:,:,2) = 0.6.*ones(size(X,1));
CO(:,:,3) = 0.6.*ones(size(X,1));
figure,surf(0.99*h*X,0.99*h*Y,0.99*h*Z,CO,'FaceAlpha',0.95),lighting phong, shading interp, hold on

h1 = J*om1'; h2 = J*om2'; h3 = J*om3';
plot3(h1(1,:),h1(2,:),h1(3,:),'b',-h1(1,:),-h1(2,:),-h1(3,:),'b','LineWidth',2)
plot3(h2(1,:),h2(2,:),h2(3,:),'k',-h2(1,:),-h2(2,:),-h2(3,:),'k','LineWidth',2)
plot3(h3(1,:),h3(2,:),h3(3,:),'r',-h3(1,:),-h3(2,:),-h3(3,:),'r','LineWidth',2)

h01 = J*om_vec1;
h02 = J*om_vec2;
h03 = J*om_vec3;
plot3(h01(1,:),h01(2,:),h01(3,:),'k*','LineWidth',2)
plot3(h02(1,:),h02(2,:),h02(3,:),'k*','LineWidth',2)
plot3(h03(1,:),h03(2,:),h03(3,:),'k*','LineWidth',2)
plot3(-h01(1,:),-h01(2,:),-h01(3,:),'k*','LineWidth',2)
plot3(-h02(1,:),-h02(2,:),-h02(3,:),'k*','LineWidth',2)
plot3(-h03(1,:),-h03(2,:),-h03(3,:),'k*','LineWidth',2) 
axis equal

%% Safe Mode Rotation (about Intermediate Axis [Body Frame z-axis])
% Superspin and Dynamic Balance
J_eff = 1.2*J33;
om0_equil = -om_vec2;
rho = (J_eff - J22)*om0_equil;
h_rho = norm(J*om0_equil + rho);

% Integrate Gyrostat Dynamics
[t, om_equil] = ode45(@(t,om)HSTspin2(t,om,J,rho), tspan, om0_equil, opts);
h_equil = J*om_equil' + rho;

om04 = -[1;10;1];
om04(2) = -sqrt((om_max*J33)^2 - (J11*om04(1))^2 - (J33*om04(3))^2)/J22;
[t, om4] = ode45(@(t,om)HSTspin2(t,om,J,rho), tspan, om04, opts);
h4 = J*om4' + rho;

[X,Y,Z] = sphere(50);
CO(:,:,1) = 0.6.*ones(size(X,1));
CO(:,:,2) = 0.6.*ones(size(X,1));
CO(:,:,3) = 0.6.*ones(size(X,1));
figure,surf(0.99*h_rho*X,0.99*h_rho*Y,0.99*h_rho*Z,CO,'FaceAlpha',0.95),lighting phong, shading interp, hold on
plot3(h_equil(1,:),h_equil(2,:),h_equil(3,:),'k*','LineWidth',2)
plot3(h4(1,:),h4(2,:),h4(3,:),'b','LineWidth',2)
axis equal

%%