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
        
        [t1,wt1] = ode45(@(t,w) HSTdynamics(t,w,mu,J,rho), tspan, w1);
        [t2,wt2] = ode45(@(t,w) HSTdynamics(t,w,mu,J,rho), tspan, w2);
        [t3,wt3] = ode45(@(t,w) HSTdynamics(t,w,mu,J,rho), tspan, w3);
        
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
[t, om4] = ode45(@(t,om)HSTdynamics(t,om,mu,J,rho), tspan, om04, opts);
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

%% Model Full HST Dynamics
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
[r_eci, v_eci] = OE2ECI(a, e, ii, RAAN, w, anom, mu);
x0 = [r_eci;v_eci];

% ODE45 Solver for Orbital Dynamics
tol = 1e-6;
opts  = odeset('reltol', tol, 'abstol', tol);
tspan = [0:10:3*T];
[t, x] = ode45(@(t,x)HSTdynamics(t,x,mu,J,rho), tspan, x0, opts);

% Plot Orbit
figure, hold on
plot3(x(:,1),x(:,2),x(:,3),'r','LineWidth',2)
axis equal
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')
title('Unperturbed HST Orbit')
earthPlot(1), axis equal
hold off

% Initial Attitude Conditions
q0 = [0;0;0;1]; % Identity Quaternion
om0 = [1;1;10];
x0 = [om0;q0];

% ODE45 Solver for Attitude Dynamics
dt = 0.01;
tspan = 0:dt:0.001*T;
[t, x] = ode45(@(t,x)HSTdynamics(t,x,mu,J,rho), tspan, x0, opts);

% Plot Quaternion
q = x(:,4:7);
figure
subplot(4,1,1),plot(t,q(:,1)),ylabel('q_1'),xlabel('t'),title('Perturbed HST Quaternion')
subplot(4,1,2),plot(t,q(:,2)),ylabel('q_2'),xlabel('t')
subplot(4,1,3),plot(t,q(:,3)),ylabel('q_3'),xlabel('t')
subplot(4,1,4),plot(t,q(:,4)),ylabel('q_4'),xlabel('t')

%% HST Static Attitude Estimation (Using Magnetometer + GPS)
stddev = 0.01;
V = stddev^2*eye(3);

rN1 = 2*rand(3,1)-1; rN1 = rN1/norm(rN1);
rN2 = 2*rand(3,1)-1; rN2 = cross(rN1,rN2); rN2 = rN2/norm(rN2);
rN = [rN1,rN2];

[y(:,:,:),~] = HSTmeasure(x,V,1,rN);
for ii = 1:length(x)
    % Generate Static Measurements  
    rB = [y(:,ii,1) y(:,ii,2)];
    
    [q_triad(:,ii),~] = triad(rN,rB);
    [q_dav(:,ii),~] = qdaven(rN,rB,stddev);
    [q_svd(:,ii),~] = qsvd(rN,rB,stddev);
    q_triad(:,ii) = sign(x(ii,4:7)').*abs(q_triad(:,ii));
    error_triad(:,ii) = norm(quat2phi(qmult(qconj(x(ii,4:7)'))*q_triad(:,ii)));
    q_dav(:,ii) = sign(x(ii,4:7)').*abs(q_dav(:,ii));
    error_dav(:,ii) = norm(quat2phi(qmult(qconj(x(ii,4:7)'))*q_dav(:,ii)));
    q_svd(:,ii) = sign(x(ii,4:7)').*abs(q_svd(:,ii));
    error_svd(:,ii) = norm(quat2phi(qmult(qconj(x(ii,4:7)'))*q_svd(:,ii)));
end
mean(error_triad),mean(error_dav),mean(error_svd)
figure
subplot(4,1,1),plot(t,q(:,1),t,q_triad(1,:),'r-'),ylabel('q_1'),xlabel('t'),title('Truth vs. Triad Quaternion')
subplot(4,1,2),plot(t,q(:,2),t,q_triad(2,:),'r-'),ylabel('q_2'),xlabel('t')
subplot(4,1,3),plot(t,q(:,3),t,q_triad(3,:),'r-'),ylabel('q_3'),xlabel('t')
subplot(4,1,4),plot(t,q(:,4),t,q_triad(4,:),'r-'),ylabel('q_4'),xlabel('t')
legend('Truth','Triad')
figure
subplot(4,1,1),plot(t,q(:,1),t,q_dav(1,:),'r-'),ylabel('q_1'),xlabel('t'),title('Truth vs. Davenport Quaternion')
subplot(4,1,2),plot(t,q(:,2),t,q_dav(2,:),'r-'),ylabel('q_2'),xlabel('t')
subplot(4,1,3),plot(t,q(:,3),t,q_dav(3,:),'r-'),ylabel('q_3'),xlabel('t')
subplot(4,1,4),plot(t,q(:,4),t,q_dav(4,:),'r-'),ylabel('q_4'),xlabel('t')
legend('Truth','Davenport')
figure
subplot(4,1,1),plot(t,q(:,1),t,q_svd(1,:),'r-'),ylabel('q_1'),xlabel('t'),title('Truth vs. SVD Quaternion')
subplot(4,1,2),plot(t,q(:,2),t,q_svd(2,:),'r-'),ylabel('q_2'),xlabel('t')
subplot(4,1,3),plot(t,q(:,3),t,q_svd(3,:),'r-'),ylabel('q_3'),xlabel('t')
subplot(4,1,4),plot(t,q(:,4),t,q_svd(4,:),'r-'),ylabel('q_4'),xlabel('t')
legend('Truth','SVD')

%% Recursive Attitude Estimation
rN = [rN1, rN2];
rB1 = y(:,:,1); rB2 = y(:,:,2);
rB = [rB1(:,1) rB2(:,1)];
q0 = triad(rN,rB);
x0 = [q0;0;0;0]; % [quaternion; gyro bias] Initialize with no bias
P0 = (10*pi/180)^2*eye(6); %10 deg. and 10 deg/sec 1-sigma uncertainty

W = 0.0001*eye(6);
V_MEKF = diag([(0.035^2)*ones(6,1); (0.00174533^2)*ones(3,1)]);
yhist = [rB1; rB2];
dt = .01;

% Obtain Noisy State with Gyro Bias
gyro_rnd = 0.1*(pi/180); % Rate Noise Density [deg/s/sqrt(Hz)]
gyro_arw = 0.1*(pi/180); % Angle Random Walk, converted from [deg/sqrt(hr)]
V_rnd = eye(3)*(gyro_rnd^2);
V_arw = eye(3)*(gyro_arw^2);
V_gyro(:,:,1) = V_rnd;
V_gyro(:,:,2) = V_arw;
[y_gyro,v_gyro] = HSTmeasure(x,V_gyro,3,[]); % Gyroscope
whist = y_gyro(1:3,:);
bias = y_gyro(4:6,:);
whist = x(:,1:3)';

% Obtain Noisy Star Tracker Measurements
V_st = diag([(0.000174533^2)*ones(1,3)]); % Arcseconds to Radians
[qhist,v_st] = HSTmeasure(x,V_st,2,[]); % Noisy Quaternion

% Run MEKF
[xhist, Phist] = mekf(x0, P0, W, V_MEKF, rN, whist, yhist, qhist, dt);

%Calculate error quaternions
e = zeros(size(q'));
for k = 1:size(q,2)
    e(:,k) = qmult(qconj(q(k,:)'))*xhist(1:4,k);
end

% Plot Attitude
figure,subplot(4,1,1),plot(t,q(:,1),t,xhist(1,:),'r-');
title('Attitude'); legend('True', 'Estimated'),xlabel('Time (s)'),ylabel('q_1')
subplot(4,1,2),plot(t,q(:,2),t,xhist(2,:),'r-'),xlabel('Time (s)'),ylabel('q_2')
subplot(4,1,3),plot(t,q(:,3),t,xhist(3,:),'r-'),xlabel('Time (s)'),ylabel('q_3')
subplot(4,1,4),plot(t,q(:,4),t,xhist(4,:),'r-'),xlabel('Time (s)'),ylabel('q_4')

% Plot Bias
figure,subplot(3,1,1),plot(t,bias(1,:),t,xhist(5,:),'r-');
title('Bias'); legend('True', 'Estimated');
subplot(3,1,2),plot(t,bias(2,:),t,xhist(6,:),'r-');
subplot(3,1,3),plot(t,bias(3,:),t,xhist(7,:),'r-');

figure;
subplot(3,1,1);
plot((180/pi)*e(1,:));
hold on
plot((180/pi)*sqrt(squeeze(Phist(1,1,:))),'r');
plot(-(180/pi)*sqrt(squeeze(Phist(1,1,:))),'r');
title('Attitude Error'); axis([0 t(end)*100 -1 1]);
subplot(3,1,2);
plot((180/pi)*e(2,:));
hold on
plot((180/pi)*sqrt(squeeze(Phist(2,2,:))),'r');
plot(-(180/pi)*sqrt(squeeze(Phist(2,2,:))),'r');
ylabel('degrees'); axis([0 t(end)*100 -1 1]);
subplot(3,1,3);
plot((180/pi)*e(3,:));
hold on
plot((180/pi)*sqrt(squeeze(Phist(3,3,:))),'r');
plot(-(180/pi)*sqrt(squeeze(Phist(3,3,:))),'r');
axis([0 t(end)*100 -1 1]);

figure;
subplot(3,1,1);
plot(xhist(5,:)-bias(1,:));
hold on
plot(2*sqrt(squeeze(Phist(4,4,:))),'r');
plot(-2*sqrt(squeeze(Phist(4,4,:))),'r');
title('Bias Error');
subplot(3,1,2);
plot(xhist(6,:)-bias(2,:));
hold on
plot(2*sqrt(squeeze(Phist(5,5,:))),'r');
plot(-2*sqrt(squeeze(Phist(5,5,:))),'r');
subplot(3,1,3);
plot(xhist(7,:)-bias(3,:));
hold on
plot(2*sqrt(squeeze(Phist(6,6,:))),'r');
plot(-2*sqrt(squeeze(Phist(6,6,:))),'r');

%% Static Estimation Monte Carlo Simulation
stddev = 0.01;
V = stddev^2*eye(3);

for jj = 1:2000
rN1 = 2*rand(3,1)-1; rN1 = rN1/norm(rN1);
rN2 = 2*rand(3,1)-1; rN2 = cross(rN1,rN2); rN2 = rN2/norm(rN2);
rN = [rN1,rN2];

[y(:,:,:),~] = HSTmeasure(x,V,1,rN);
for ii = 1:length(x)
    % Generate Static Measurements  
    rB = [y(:,ii,1) y(:,ii,2)];
    
    [q_triad(:,ii),~] = triad(rN,rB);
    [q_dav(:,ii),~] = qdaven(rN,rB,stddev);
    [q_svd(:,ii),~] = qsvd(rN,rB,stddev);
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
end
figure
subplot(3,1,1),plot(1:1:jj,triad_mean,'.'),xlabel('N'),ylabel('Error (deg)')
title('TRIAD Algorithm Error')
subplot(3,1,2),plot(1:1:jj,dav_mean,'.'),xlabel('N'),ylabel('Error (deg)')
title('Davenport q-method Error')
subplot(3,1,3),plot(1:1:jj,svd_mean,'.'),xlabel('N'),ylabel('Error (deg)')
title('SVD Algorithm Error')