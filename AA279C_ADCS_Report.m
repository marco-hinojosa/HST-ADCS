%% AA279C ADCS Report
% Marco Hinojosa
% 06181747
close all,clear all,clc

%% Superspin and Dynamic Balance (All Calcs Done in Principal Frame)
% Principal Axes of Inertia
J11 = 28525.53;
J22 = 174815.86;
J33 = 181630.81;
J = diag([J11 J22 J33]); % [x-axis z-axis y-axis]

% Initialize Parameters and Nominal Spin Conditions
om_max = 10; % RPM
h = J33*om_max;
rho = 0;
om_vec1 = [om_max;0;0];
om_vec1(1) = sqrt((h)^2 - (J22*om_vec1(2))^2 - (J33*om_vec1(3))^2)/J11;
om_vec2 = [0;om_max;0];
om_vec2(2) = sqrt((h)^2 - (J11*om_vec2(1))^2 - (J33*om_vec2(3))^2)/J22;
om_vec3 = [0;0;om_max];

% Compute Perturbed Attitude Dynamics
tol = 1e-6; opts  = odeset('reltol', tol, 'abstol', tol);
tspan = [0:0.01:5];

count = 1; c = 4;
for ii = 0:c
    for jj = 0:c       
        w1 = [0;ii;jj];
        w1(1) = sqrt((h)^2-(J22*w1(2))^2-(J33*w1(3))^2)/J11;
        w2 = [ii;0;jj];
        w2(2) = sqrt((h)^2-(J11*w2(1))^2-(J33*w2(3))^2)/J22;
        w3 = [0;ii;jj];
        w3(3) = sqrt((h)^2-(J22*w3(2))^2-(J11*w3(1))^2)/J33;
        
        [t1,wt1] = ode45(@(t,w) HSTspin(t,w,J,rho), tspan, w1);
        [t2,wt2] = ode45(@(t,w) HSTspin(t,w,J,rho), tspan, w2);
        [t3,wt3] = ode45(@(t,w) HSTspin(t,w,J,rho), tspan, w3);
        
        ht1(:,:,count) = J*wt1';
        ht2(:,:,count) = J*wt2';
        ht3(:,:,count) = J*wt3';
        
        wt1(:,:,count) = wt1;
        wt2(:,:,count) = wt2;
        wt3(:,:,count) = wt3;
        
        count = count + 1;
    end
end

% Plot Perturbed Dynamics
figure,hold on
for ii = 1:count-1
    plot3(ht1(1,:,ii),ht1(2,:,ii),ht1(3,:,ii),'r','LineWidth',1)
    plot3(ht2(1,:,ii),ht2(2,:,ii),ht2(3,:,ii),'b','LineWidth',1)
    plot3(ht3(1,:,ii),ht3(2,:,ii),ht3(3,:,ii),'k','LineWidth',1)
    plot3(-ht1(1,:,ii),-ht1(2,:,ii),-ht1(3,:,ii),'r','LineWidth',1)
    plot3(-ht2(1,:,ii),-ht2(2,:,ii),-ht2(3,:,ii),'b','LineWidth',1)
    plot3(-ht3(1,:,ii),-ht3(2,:,ii),-ht3(3,:,ii),'k','LineWidth',1)
end

% Plot Momentum Sphere
[X,Y,Z] = sphere(50);
CO(:,:,1) = 0.6.*ones(size(X,1));
CO(:,:,2) = 0.6.*ones(size(X,1));
CO(:,:,3) = 0.6.*ones(size(X,1));
surf(0.99*h*X,0.99*h*Y,0.99*h*Z,CO,'FaceAlpha',0.95),lighting phong, shading interp

% Plot Equilibrium Points
hequil_1 = J*om_vec1;
hequil_2 = J*om_vec2;
hequil_3 = J*om_vec3;
plot3(hequil_1(1,:),hequil_1(2,:),hequil_1(3,:),'k*','LineWidth',2)
plot3(hequil_2(1,:),hequil_2(2,:),hequil_2(3,:),'k*','LineWidth',2)
plot3(hequil_3(1,:),hequil_3(2,:),hequil_3(3,:),'k*','LineWidth',2)
plot3(-hequil_1(1,:),-hequil_1(2,:),-hequil_1(3,:),'k*','LineWidth',2)
plot3(-hequil_2(1,:),-hequil_2(2,:),-hequil_2(3,:),'k*','LineWidth',2)
plot3(-hequil_3(1,:),-hequil_3(2,:),-hequil_3(3,:),'k*','LineWidth',2)
title('Momentum Sphere w/o Gyrostat')
legend('Body X-Axis (Major)','Body Z-Axis (Inter)','Body Y-Axis (Minor)')
axis equal

% ---------------------------------------------------------------------
% Safe Mode Rotation (about Intermediate Axis [Body Frame Y-Axis])
J_eff = 1.2*J33;
om0_equil = om_vec3;
rho = (J_eff - J33)*om0_equil;
h_rho = norm(J*om0_equil + rho);

% Integrate Gyrostat Dynamics
[t, om_equil] = ode45(@(t,om)HSTspin(t,om,J,rho), tspan, om0_equil, opts);
h_equil = J*om_equil' + rho;

om04 = [2;2;10];
om04(3) = sqrt((om_max*J33)^2 - (J11*om04(1))^2 - (J22*om04(2))^2)/J33;
[t, om4] = ode45(@(t,om)HSTspin(t,om,J,rho), tspan, om04, opts);
h4 = J*om4' + rho;

[X,Y,Z] = sphere(50);
CO(:,:,1) = 0.6.*ones(size(X,1));
CO(:,:,2) = 0.6.*ones(size(X,1));
CO(:,:,3) = 0.6.*ones(size(X,1));
figure,hold on

% Plot Perturbed Dynamics
plot3(h4(1,:),h4(2,:),h4(3,:),'b','LineWidth',2)

% Plot Momentum Sphere
surf(0.99*h_rho*X,0.99*h_rho*Y,0.99*h_rho*Z,CO,'FaceAlpha',0.95),lighting phong, shading interp

% Plot Equilibrium Point
plot3(h_equil(1,:),h_equil(2,:),h_equil(3,:),'k*','LineWidth',2)
title('Momentum Sphere w/ Gyrostat')
legend('Stable Spin about Minor Axis')
axis equal

% %% Safe Mode Rotation (about Intermediate Axis [Body Frame Negative Z-Axis])
% J_eff = 1.2*J33;
% om0_equil = -om_vec2;
% rho = (J_eff - J22)*om0_equil;
% h_rho = norm(J*om0_equil + rho);
% 
% % Integrate Gyrostat Dynamics
% [t, om_equil] = ode45(@(t,om)HSTspin(t,om,J,rho), tspan, om0_equil, opts);
% h_equil = J*om_equil' + rho;
% 
% om04 = -[2;10;2];
% om04(2) = -sqrt((om_max*J33)^2 - (J11*om04(1))^2 - (J33*om04(3))^2)/J22;
% [t, om4] = ode45(@(t,om)HSTspin(t,om,J,rho), tspan, om04, opts);
% h4 = J*om4' + rho;
% 
% [X,Y,Z] = sphere(50);
% CO(:,:,1) = 0.6.*ones(size(X,1));
% CO(:,:,2) = 0.6.*ones(size(X,1));
% CO(:,:,3) = 0.6.*ones(size(X,1));
% figure,plot3(h4(1,:),h4(2,:),h4(3,:),'b','LineWidth',2),hold on
% plot3(h_equil(1,:),h_equil(2,:),h_equil(3,:),'k*','LineWidth',2)
% surf(0.99*h_rho*X,0.99*h_rho*Y,0.99*h_rho*Z,CO,'FaceAlpha',0.95),lighting phong, shading interp
% title('Momentum Sphere w/ Gyrostat')
% legend('Body Negative Z-Axis (Inter)')
% axis equal

%% Model HST Orbital and Attitude Dynamics
% Reference Orbital Elements and Parameters
a = 6917.5; % Orbit Semi-Major Axis (km)
e = 0.000287; % Eccentricity
ii = 28.47; % Inclination (deg)
RAAN = 176.23; % Right Ascension of Ascending Node (deg)
w = 82.61; % Argument of Perigee (deg)
anom = 319.41; % Mean Anomaly (deg)
mu = 3.986e5; % Earth Standard Gravitational Parameter (km^3/s^2)
T = 2*pi*sqrt((a^3)/mu); % Orbital Period (s)

% Convert to Earth-Centered Inertial Frame Coordinates
[r_eci, v_eci] = OE2ECI(a, e, ii, RAAN, w, anom, mu);
x0 = [r_eci;v_eci;0;0;0;0;0;0;0];

% ODE45 Solver for Orbital Dynamics
tol = 1e-6;
opts  = odeset('reltol', tol, 'abstol', tol);
tspan = [0:10:3*T];
[t, x] = ode45(@(t,x)HSTdynamics(mu,J,rho,x), tspan, x0, opts);

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

% Initial Attitude Conditions
Q0 = [0 -1 0;1 0 0;0 0 1];
q0 = Q2quat(Q0);
om0 = [1;1;10];
x0 = [r_eci;v_eci;om0;q0];

% ODE45 Solver for Attitude Dynamics
tspan = 0:0.001:0.01*T;
[t, x] = ode45(@(t,x)HSTdynamics(mu,J,rho,x), tspan, x0, opts);

% Plot Quaternion
q = x(:,10:12);
figure,plot3(q(:,1),q(:,2),q(:,3),'LineWidth',0.75)
title('Perturbed HST Quaternion')
