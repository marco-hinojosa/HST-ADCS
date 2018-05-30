%% AA279C ADCS Report
% Marco Hinojosa
% 06181747
close all,clear all,clc

%% Superspin and Dynamic Balance (All Calcs Done in Principal Frame)
% Reference Orbital Elements and Parameters
a = 6917.5; % Orbit Semi-Major Axis (km)
e = 0.000287; % Eccentricity
i = 28.47; % Inclination (deg)
RAAN = 176.23; % Right Ascension of Ascending Node (deg)
w = 82.61; % Argument of Perigee (deg)
anom = 319.41; % Mean Anomaly (deg)
mu = 3.986e5; % Earth Standard Gravitational Parameter (km^3/s^2)
T = 2*pi*sqrt((a^3)/mu); % Orbital Period (s)

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
        
        [t1,wt1] = ode45(@(t,w) HSTdynamics(t,w,mu,J,rho,0), tspan, w1);
        [t2,wt2] = ode45(@(t,w) HSTdynamics(t,w,mu,J,rho,0), tspan, w2);
        [t3,wt3] = ode45(@(t,w) HSTdynamics(t,w,mu,J,rho,0), tspan, w3);
        
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

% -------------------------------------------------------------------------
% Safe Mode Rotation (about Intermediate Axis [Body Frame Y-Axis])
J_eff = 1.2*J33;
om0_equil = om_vec3;
rho = (J_eff - J33)*om0_equil;
h_rho = norm(J*om0_equil + rho);

% Integrate Gyrostat Dynamics
[t, om_equil] = ode45(@(t,om)HSTdynamics(t,om,mu,J,rho), tspan, om0_equil, opts);
h_equil = J*om_equil' + rho*ones(1,length(om_equil));

om04 = [2;2;10];
om04(3) = sqrt((om_max*J33)^2 - (J11*om04(1))^2 - (J22*om04(2))^2)/J33;
[t, om4] = ode45(@(t,om)HSTdynamics(t,om,mu,J,rho,0), tspan, om04, opts);
h4 = J*om4' + rho*ones(1,length(om4));

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

% ODE45 Solver for Orbital Dynamics
tol = 1e-6;
opts  = odeset('reltol', tol, 'abstol', tol);
tspan = [0:10:3*T];
[t, x] = ode45(@(t,x)HSTdynamics(t,x,mu,J,rho,0), tspan, x0, opts);

% Plot Orbit
figure, hold on
plot3(x(:,1),x(:,2),x(:,3),'r','LineWidth',2)
axis equal
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')
title('HST Orbit Dynamics')
earthPlot(1), axis equal
hold off

%% Model HST Attitude Dynamics
clc,clearvars
mu = 3.986e5; % Earth Standard Gravitational Parameter (km^3/s^2)
% Principal Axes of Inertia
J11 = 28525.53;
J22 = 174815.86;
J33 = 181630.81;
J = diag([J11 J22 J33]); % [x-axis z-axis y-axis]
rho = 1.0e+05*[0;0;3.6326]; % From Superspin

% Initial Attitude Conditions
q0 = [0;0;0;1]; % Identity Quaternion
om0 = [0.25;0.25;2.5]; % rad/s
x0 = [om0;q0];

% ODE45 Solver for Attitude Dynamics
dt = 0.01;
tspan = 0:dt:60;
tol = 1e-6; opts  = odeset('reltol', tol, 'abstol', tol);
[t, x] = ode45(@(t,x)HSTdynamics(t,x,mu,J,rho,0), tspan, x0, opts);

% Plot Quaternion
om_true = x(:,1:3); q_true = x(:,4:7);
figure
subplot(4,1,1),plot(t,q_true(:,1)),ylabel('q_1'),xlabel('t'),title('HST Quaternion Dynamics')
subplot(4,1,2),plot(t,q_true(:,2)),ylabel('q_2'),xlabel('t')
subplot(4,1,3),plot(t,q_true(:,3)),ylabel('q_3'),xlabel('t')
subplot(4,1,4),plot(t,q_true(:,4)),ylabel('q_4'),xlabel('t')

%% Sensor Specifications
% Star Tracker Specifications (Ball Aerospace CT-601)
error_ST = 300; % Arcseconds
error_ST = error_ST*(1/3600)*(pi/180); % Radians
V_st = diag((error_ST^2)*ones(1,3));
%V_st = diag(0.0003*ones(1,3)); % In-Class Example

% Gyro Specifications (HST Gas-Bearing Gyros)
gyro_rnd = 0.1*(pi/180); % Rate Noise Density, converted from [deg/s/sqrt(Hz)]
gyro_arw = 0.1*(pi/180); % Angle Random Walk, converted from [deg/sqrt(hr)]
W_rnd = (gyro_rnd^2)*eye(3); % Gyro Covariance
W_arw = (gyro_arw^2)*eye(3); % Bias Covariance
%W_rnd = (0.3046)*eye(3); % In-Class Example
%W_arw = (0.3046)*eye(3); % In-Class Example
W_gyro = blkdiag(W_rnd, W_arw);

% Magnetometer Specifications
error_mag = 4*(pi/180); % Radians
V_mag = (error_mag^2)*eye(3);
%V_mag = (0.0076)*eye(3); % In-Class Example

%% Static Attitude Estimation (Using Magnetometer + GPS)
% Noisy Magnetometer Body-Frame Measurements
rN1 = 2*rand(3,1)-1; rN1 = rN1/norm(rN1);
rN2 = 2*rand(3,1)-1; rN2 = cross(rN1,rN2); rN2 = rN2/norm(rN2);
rN = [rN1,rN2];
[rB_noisy,~] = HSTmeasure(x,V_mag,1,rN);

% Perform Static Estimation at Each Time Step
for ii = 1:length(x)
    rB = [rB_noisy(1:3,ii) rB_noisy(4:6,ii)];   
    [q_triad(:,ii),~] = triad(rN,rB);
    [q_dav(:,ii),~] = qdaven(rN,rB,error_mag);
    [q_svd(:,ii),~] = qsvd(rN,rB,error_mag);
    q_triad(:,ii) = sign(x(ii,4:7)').*abs(q_triad(:,ii));
    error_triad(:,ii) = norm(quat2phi(qmult(qconj(x(ii,4:7)'))*q_triad(:,ii)));
    q_dav(:,ii) = sign(x(ii,4:7)').*abs(q_dav(:,ii));
    error_dav(:,ii) = norm(quat2phi(qmult(qconj(x(ii,4:7)'))*q_dav(:,ii)));
    q_svd(:,ii) = sign(x(ii,4:7)').*abs(q_svd(:,ii));
    error_svd(:,ii) = norm(quat2phi(qmult(qconj(x(ii,4:7)'))*q_svd(:,ii)));
end

% Compute Mean Error and Plot Results
mean(error_triad)*180/pi,mean(error_dav)*180/pi,mean(error_svd)*180/pi
figure
subplot(4,1,1),plot(t,q_true(:,1),t,q_triad(1,:),'r-'),ylabel('q_1'),xlabel('t'),title('Truth vs. Triad Quaternion')
subplot(4,1,2),plot(t,q_true(:,2),t,q_triad(2,:),'r-'),ylabel('q_2'),xlabel('t')
subplot(4,1,3),plot(t,q_true(:,3),t,q_triad(3,:),'r-'),ylabel('q_3'),xlabel('t')
subplot(4,1,4),plot(t,q_true(:,4),t,q_triad(4,:),'r-'),ylabel('q_4'),xlabel('t')
legend('Truth','Triad')
figure
subplot(4,1,1),plot(t,q_true(:,1),t,q_dav(1,:),'r-'),ylabel('q_1'),xlabel('t'),title('Truth vs. Davenport Quaternion')
subplot(4,1,2),plot(t,q_true(:,2),t,q_dav(2,:),'r-'),ylabel('q_2'),xlabel('t')
subplot(4,1,3),plot(t,q_true(:,3),t,q_dav(3,:),'r-'),ylabel('q_3'),xlabel('t')
subplot(4,1,4),plot(t,q_true(:,4),t,q_dav(4,:),'r-'),ylabel('q_4'),xlabel('t')
legend('Truth','Davenport')
figure
subplot(4,1,1),plot(t,q_true(:,1),t,q_svd(1,:),'r-'),ylabel('q_1'),xlabel('t'),title('Truth vs. SVD Quaternion')
subplot(4,1,2),plot(t,q_true(:,2),t,q_svd(2,:),'r-'),ylabel('q_2'),xlabel('t')
subplot(4,1,3),plot(t,q_true(:,3),t,q_svd(3,:),'r-'),ylabel('q_3'),xlabel('t')
subplot(4,1,4),plot(t,q_true(:,4),t,q_svd(4,:),'r-'),ylabel('q_4'),xlabel('t')
legend('Truth','SVD')

%% Recursive Attitude Estimation
% Noisy Magnetometer Body-Frame Measurements
rN1 = 2*rand(3,1)-1; rN1 = rN1/norm(rN1);
rN2 = 2*rand(3,1)-1; rN2 = cross(rN1,rN2); rN2 = rN2/norm(rN2);
rN = [rN1,rN2];
[rB_noisy,~] = HSTmeasure(x,V_mag,1,rN);

% Noisy Star Tracker Quaternion Measurements
[q_noisy,~] = HSTmeasure(x,V_st,2,[]);

% Noisy Gyro State Measurements (with Bias)
[om_noisy,bias] = HSTmeasure(x,W_gyro,3,[]);

% MEKF Noisy Measurement Histories
whist = om_noisy;
yhist = [q_noisy; rB_noisy];

% MEKF Initial Values
q0 = q_true(1,:)';
x0 = [q0;0;0;0]; % [quaternion; gyro bias]
P0 = (10*pi/180)^2*eye(6);
W_MEKF = W_gyro;
V_MEKF = blkdiag(V_st, V_mag, V_mag);

% Run MEKF
[xhist, Phist] = mekf(x0, P0, W_MEKF, V_MEKF, rN, whist, yhist, dt);

% Calculate Error Quaternions
e = zeros(3,size(q_true,1));
for k = 1:size(q_true,1)
    e(:,k) = quat2phi(qmult(qconj(q_true(k,:)'))*xhist(1:4,k));
end

% Plot Attitude
figure,
subplot(4,1,1),plot(t,q_true(:,1),t,xhist(1,:),'r-'),xlabel('Time (s)'),ylabel('q_1')
title('Attitude'); legend('True', 'Estimated')
subplot(4,1,2),plot(t,q_true(:,2),t,xhist(2,:),'r-'),xlabel('Time (s)'),ylabel('q_2')
subplot(4,1,3),plot(t,q_true(:,3),t,xhist(3,:),'r-'),xlabel('Time (s)'),ylabel('q_3')
subplot(4,1,4),plot(t,q_true(:,4),t,xhist(4,:),'r-'),xlabel('Time (s)'),ylabel('q_4')

% Plot Bias
figure,
subplot(3,1,1),plot(t,bias(1,:),t,xhist(5,:),'r-');
title('Bias'); legend('True', 'Estimated');
subplot(3,1,2),plot(t,bias(2,:),t,xhist(6,:),'r-');
subplot(3,1,3),plot(t,bias(3,:),t,xhist(7,:),'r-');
xlabel('Time (s)')

% Plot Errors
figure;
subplot(3,1,1);
plot(t,(180/pi)*e(1,:)); axis([0 t(end) -0.25 0.25]); hold on
plot(t,2*(180/pi)*sqrt(squeeze(Phist(1,1,:))),'r');
plot(t,-2*(180/pi)*sqrt(squeeze(Phist(1,1,:))),'r');
title('Attitude Error');
subplot(3,1,2);
plot(t,(180/pi)*e(2,:)); axis([0 t(end) -0.25 0.25]); hold on
plot(t,2*(180/pi)*sqrt(squeeze(Phist(2,2,:))),'r');
plot(t,-2*(180/pi)*sqrt(squeeze(Phist(2,2,:))),'r');
ylabel('Degrees');
subplot(3,1,3);
plot(t,(180/pi)*e(3,:)); axis([0 t(end) -0.25 0.25]); hold on
plot(t,2*(180/pi)*sqrt(squeeze(Phist(3,3,:))),'r');
plot(t,-2*(180/pi)*sqrt(squeeze(Phist(3,3,:))),'r');
xlabel('Time (s)')

figure;
subplot(3,1,1);
plot(t,xhist(5,:)-bias(1,:)); axis([0 t(end) -0.25 0.25]); hold on
plot(t,2*sqrt(squeeze(Phist(4,4,:))),'r');
plot(t,-2*sqrt(squeeze(Phist(4,4,:))),'r');
title('Bias Error');
subplot(3,1,2);
plot(t,xhist(6,:)-bias(2,:)); axis([0 t(end) -0.25 0.25]); hold on
plot(t,2*sqrt(squeeze(Phist(5,5,:))),'r');
plot(t,-2*sqrt(squeeze(Phist(5,5,:))),'r');
subplot(3,1,3);
plot(t,xhist(7,:)-bias(3,:)); axis([0 t(end) -0.25 0.25]); hold on
plot(t,2*sqrt(squeeze(Phist(6,6,:))),'r');
plot(t,-2*sqrt(squeeze(Phist(6,6,:))),'r');
xlabel('Time (s)')

%% Static Estimation Monte Carlo Simulation
for jj = 1:2000
    % Noisy Magnetometer Body-Frame Measurements
    rN1 = 2*rand(3,1)-1; rN1 = rN1/norm(rN1);
    rN2 = 2*rand(3,1)-1; rN2 = cross(rN1,rN2); rN2 = rN2/norm(rN2);
    rN = [rN1,rN2];
    [rB_noisy,~] = HSTmeasure(x,V_mag,1,rN);
    for ii = 1:length(x)
        rB = [rB_noisy(1:3,ii) rB_noisy(4:6,ii)];
        [q_triad(:,ii),~] = triad(rN,rB);
        [q_dav(:,ii),~] = qdaven(rN,rB,error_mag);
        [q_svd(:,ii),~] = qsvd(rN,rB,error_mag);
        q_triad(:,ii) = sign(x(ii,4:7)').*abs(q_triad(:,ii));
        error_triad(:,ii) = norm(quat2phi(qmult(qconj(x(ii,4:7)'))*q_triad(:,ii)));
        q_dav(:,ii) = sign(x(ii,4:7)').*abs(q_dav(:,ii));
        error_dav(:,ii) = norm(quat2phi(qmult(qconj(x(ii,4:7)'))*q_dav(:,ii)));
        q_svd(:,ii) = sign(x(ii,4:7)').*abs(q_svd(:,ii));
        error_svd(:,ii) = norm(quat2phi(qmult(qconj(x(ii,4:7)'))*q_svd(:,ii)));
    end
    triad_mean(jj) = mean(error_triad)*180/pi;
    dav_mean(jj) = mean(error_dav)*180/pi;
    svd_mean(jj) = mean(error_svd)*180/pi;
    jj
end

% Plot Mean Error Results
figure
subplot(3,1,1),plot(1:1:jj,triad_mean,'.'),xlabel('N'),ylabel('Error (deg)')
title('TRIAD Algorithm Error')
subplot(3,1,2),plot(1:1:jj,dav_mean,'.'),xlabel('N'),ylabel('Error (deg)')
title('Davenport q-method Error')
subplot(3,1,3),plot(1:1:jj,svd_mean,'.'),xlabel('N'),ylabel('Error (deg)')
title('SVD Algorithm Error')

%% Model HST Dynamics with Environmental Disturbances
clc,clearvars
% Reference Orbital Elements and Parameters
a = 6917.5; % Orbit Semi-Major Axis (km)
e = 0.000287; % Eccentricity
i = 28.47; % Inclination (deg)
RAAN = 176.23; % Right Ascension of Ascending Node (deg)
w = 82.61; % Argument of Perigee (deg)
anom = 319.41; % Mean Anomaly (deg)
mu = 3.986e5; % Earth Standard Gravitational Parameter (km^3/s^2)
T = 2*pi*sqrt((a^3)/mu); % Orbital Period (s)

% Principal Axes of Inertia
J11 = 28525.53;
J22 = 174815.86;
J33 = 181630.81;
J = diag([J11 J22 J33]); % [x-axis z-axis y-axis]
rho = 1.0e+05*[0;0;3.6326]; % From Superspin

% Compute Environmental Disturbance Parameters (HST-Specific Values)
CG = [-0.01; 6.03; 0.00];
n1 = [1 0 0]'; n2 = [-1 0 0]'; n3 = [0 1 0]'; n4 = [0 -1 0]'; n5 = [0 0 1]'; n6 = [0 0 -1]';
c1 = [0.64; 6.59; 0.00]; c2 = [-1.00; 6.59; 0.00]; c3 = [-0.03; 10.28; 0.00];
c4 = [-0.03; 0.72; 0.00]; c5 = [-0.01; 6.28; 1.85]; c6 = [-0.01; 6.28; -1.85];

S = [122.73, 122.73, 20.97, 20.97, 54.04, 54.04]; % Surface Area of Each Side, m^2
n = [n1 n2 n3 n4 n5 n6]; % Surface Normal Vectors
c = [c1 c2 c3 c4 c5 c6] - CG*ones(1,6); % Surface Centroid Coordinates w.r.t. CG, m
params = [S;n;c];

% Initial Conditions
Q0 = [0 1 0; -1 0 0; 0 0 1];
q0 = Q2quat(Q0);
[r_eci, v_eci] = OE2ECI(a, e, i, RAAN, w, anom, mu);
%q0 = [0;0;0;1]; % Identity Quaternion
om0 = [0;0;0]; % rad/s, Start From Rest
x0 = [om0;q0;r_eci;v_eci];

% ODE45 Solver for Full HST Dynamics
dt = 3;
tspan = 0:dt:15*T;
tol = 1e-6; opts  = odeset('reltol', tol, 'abstol', tol);
[t, x] = ode45(@(t,x)HSTdynamics(t,x,mu,J,rho,params,1), tspan, x0, opts);

% Plot Orbit
figure, hold on
plot3(x(:,8),x(:,9),x(:,10),'r','LineWidth',2),axis equal
xlabel('[km]'),ylabel('[km]'),zlabel('[km]')
title('HST Orbit Dynamics')
earthPlot(1), axis equal
hold off

% Plot Quaternion
om_true = x(:,1:3);
figure
subplot(3,1,1),plot(t,om_true(:,1)),ylabel('\omega_x'),xlabel('t'),title('HST Angular Velocity with Disturbances')
subplot(3,1,2),plot(t,om_true(:,2)),ylabel('\omega_y'),xlabel('t')
subplot(3,1,3),plot(t,om_true(:,3)),ylabel('\omega_z'),xlabel('t')

% Plot Quaternion
q_true = x(:,4:7);
figure
subplot(4,1,1),plot(t,q_true(:,1)),ylabel('q_1'),xlabel('t'),title('HST Quaternion Dynamics with Disturbances')
subplot(4,1,2),plot(t,q_true(:,2)),ylabel('q_2'),xlabel('t')
subplot(4,1,3),plot(t,q_true(:,3)),ylabel('q_3'),xlabel('t')
subplot(4,1,4),plot(t,q_true(:,4)),ylabel('q_4'),xlabel('t')

% Plot Orbital Decay
orb_alt = vecnorm(x(:,8:10)');
figure,plot(t,orb_alt),title('Orbital Decay Over 15 Periods')
xlabel('Time (s)'),ylabel('Radial Distance from Earth (km)')

% Plot Angular Momentum
for ii = 1:size(x,1)
    h(ii) = norm(J*x(ii,1:3)');
end
figure;
plot(t,h,'b','LineWidth',1)
title('Angular Momentum')

%% Recover Air Drag and Torque History
for ii = 1:size(x,1)
    Q = quat2Q(x(ii,4:7)');
    [Drag(:,ii),Torque(:,ii)] = HSTdrag(x(ii,8:10)',x(ii,11:13)',Q,S,n,c);
    Tau_g(:,ii) = cross(3*mu/(dot(x(ii,8:10)',x(ii,8:10)')^(5/2))*x(ii,8:10)', Q*J/(1000^2)*x(ii,8:10)'); % [N-km]
    Tau_g(:,ii) = Q'*Tau_g(:,ii)*1000; % Rotate into Body Frame and Convert Units, [N-m]
end

max_drag = max(vecnorm(Drag));
max_dragtorque = max(vecnorm(Torque));
max_gravtorque = max(vecnorm(Tau_g));
average_torque = mean(vecnorm(Torque + Tau_g));