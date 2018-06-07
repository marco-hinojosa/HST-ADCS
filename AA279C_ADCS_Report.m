%% AA279C ADCS Report
% Marco Hinojosa
% 06181747
close all,clear all,clc

%% Initialize and Save Constant Parameters
% Orbital Elements and Parameters
a = 6917.5; % Orbit Semi-Major Axis (km)
e = 0.000287; % Eccentricity
i = 28.47; % Inclination (deg)
RAAN = 176.23; % Right Ascension of Ascending Node (deg)
w = 82.61; % Argument of Perigee (deg)
anom = 319.41; % Mean Anomaly (deg)
mu = 3.986e5; % Earth Standard Gravitational Parameter (km^3/s^2)
T = 2*pi*sqrt((a^3)/mu); % Orbital Period (s)
save('Parameters\orbital','a','e','i','RAAN','w','anom','mu','T')

% Mass Properties, Principle Axes of Inertia
J11 = 28525.53;
J22 = 174815.86;
J33 = 181630.81;
J = diag([J11 J22 J33]); % [x-axis z-axis y-axis]
save('Parameters\mass','J11','J22','J33','J')

% Sensor Specifications
% Star Tracker (Ball Aerospace CT-601)
error_ST = 300; % Arcseconds
error_ST = error_ST*(1/3600)*(pi/180); % Radians
V_st = diag((error_ST^2)*ones(1,3));
% Gyro (HST Gas-Bearing Gyros)
gyro_rnd = 0.1*(pi/180); % Rate Noise Density, converted from [deg/s/sqrt(Hz)]
gyro_arw = 0.1*(pi/180); % Angle Random Walk, converted from [deg/sqrt(hr)]
W_rnd = (gyro_rnd^2)*eye(3); % Gyro Covariance
W_arw = (gyro_arw^2)*eye(3); % Bias Covariance
W_gyro = blkdiag(W_rnd, W_arw);
% Magnetometer 
error_mag = 4*(pi/180); % Radians
V_mag = (error_mag^2)*eye(3);
save('Parameters\sensor','error_ST','V_st','gyro_rnd','gyro_arw','W_gyro','error_mag','V_mag')

% Compute Environmental Disturbance Parameters (HST-Specific Values)
CG = [-0.01; 6.03; 0.00];
n1 = [1 0 0]'; n2 = [-1 0 0]'; n3 = [0 1 0]'; n4 = [0 -1 0]'; n5 = [0 0 1]'; n6 = [0 0 -1]';
c1 = [0.64; 6.59; 0.00]; c2 = [-1.00; 6.59; 0.00]; c3 = [-0.03; 10.28; 0.00];
c4 = [-0.03; 0.72; 0.00]; c5 = [-0.01; 6.28; 1.85]; c6 = [-0.01; 6.28; -1.85];
S = [122.73, 122.73, 20.97, 20.97, 54.04, 54.04]; % Surface Area of Each Side, m^2
n = [n1 n2 n3 n4 n5 n6]; % Surface Normal Vectors
c = [c1 c2 c3 c4 c5 c6] - CG*ones(1,6); % Surface Centroid Coordinates w.r.t. CG, m
params = [S;n;c];
% HST Reaction Wheel Actuator Jacobian
Bw = [-1/sqrt(3),  1/sqrt(3),  1/sqrt(3), -1/sqrt(3);
       1/sqrt(3),  1/sqrt(3),  1/sqrt(3),  1/sqrt(3);
      -1/sqrt(3), -1/sqrt(3),  1/sqrt(3),  1/sqrt(3)];
save('Parameters\geometry','params','Bw')

%% Superspin and Dynamic Balance (All Calcs Done in Principal Frame)
clc,clearvars,close all
load('Parameters\orbital')
load('Parameters\mass')

% Initialize Nominal Spin Conditions
om_max = 10; % RPM
dt = J33*om_max;
rho = 0;
om_vec1 = [om_max;0;0];
om_vec1(1) = sqrt((dt)^2 - (J22*om_vec1(2))^2 - (J33*om_vec1(3))^2)/J11;
om_vec2 = [0;om_max;0];
om_vec2(2) = sqrt((dt)^2 - (J11*om_vec2(1))^2 - (J33*om_vec2(3))^2)/J22;
om_vec3 = [0;0;om_max];

% Compute Perturbed Attitude Dynamics
tol = 1e-6; opts  = odeset('reltol', tol, 'abstol', tol);
tspan = [0:0.01:5];

count = 1; c = 4;
for ii = 0:c
    for jj = 0:c       
        w1 = [0;ii;jj];
        w1(1) = sqrt((dt)^2-(J22*w1(2))^2-(J33*w1(3))^2)/J11;
        w2 = [ii;0;jj];
        w2(2) = sqrt((dt)^2-(J11*w2(1))^2-(J33*w2(3))^2)/J22;
        w3 = [0;ii;jj];
        w3(3) = sqrt((dt)^2-(J22*w3(2))^2-(J11*w3(1))^2)/J33;
        
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
surf(0.99*dt*X,0.99*dt*Y,0.99*dt*Z,CO,'FaceAlpha',0.95),lighting phong, shading interp

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
clc,clearvars,close all
load('Parameters\orbital')
load('Parameters\mass')

% Convert to Earth-Centered Inertial Frame Coordinates
[r_eci, v_eci] = OE2ECI(a, e, i, RAAN, w, anom, mu);
x0 = [r_eci;v_eci];
rho = 0;

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
clc,clearvars,close all
load('Parameters\orbital')
load('Parameters\mass')

% From Superspin
rho = 1.0e+05*[0;0;3.6326];

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

%% Static Attitude Estimation (Using Magnetometer + GPS)
load('Parameters\sensor')

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
triad_err = mean(error_triad)*180/pi
dav_err = mean(error_dav)*180/pi
svd_err = mean(error_svd)*180/pi
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

%% Recursive Attitude Estimation
load('Parameters\orbital')
load('Parameters\mass')

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

%% Model HST Dynamics with Environmental Disturbances
clc,clearvars,close all
load('Parameters\orbital')
load('Parameters\mass')
load('Parameters\geometry')

% From Superspin
rho = 1.0e+05*[0;0;3.6326]; 

% Initial Conditions
[r_eci, v_eci] = OE2ECI(a, e, i, RAAN, w, anom, mu);
q0 = [0;0;0;1]; % Identity Quaternion
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

% Plot Angular Velocity
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
    dt(ii) = norm(J*x(ii,1:3)');
end
figure;
plot(t,dt,'b','LineWidth',1)
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

%% Simulate LQR-Controlled Attitude Dynamics with Disturbances and MEKF
clc,clearvars,close all
load('Parameters\orbital')
load('Parameters\mass')
load('Parameters\sensor')
load('Parameters\geometry')

% Initial Conditions
[r_eci, v_eci] = OE2ECI(a, e, i, RAAN, w, anom, mu);
qstar = [1/sqrt(2); 1/sqrt(2); 0; 0]; % Desired Quaternion
phi_rand = rand(3,1);
%phi0 = (pi/2)*phi_rand/norm(phi_rand); % Normalize and Set 90° Initial Error
%phi0 = (pi/4)*phi_rand/norm(phi_rand); % Normalize and Set 45° Initial Error
%phi0 = (pi/8)*phi_rand/norm(phi_rand); % Normalize and Set 22.5° Initial Error
phi0 = (pi/16)*phi_rand/norm(phi_rand); % Normalize and Set 11.25° Initial Error
q_err = phi2quat(phi0); q(:,1) = qmult(q_err)*qstar;
om0 = [0; 0; 0]; % rad/s, Zero Initial Angular Velocity
x0 = [phi0; om0; r_eci; v_eci];

% Compute LQR Gains about Linearized Dynamics
A = [zeros(3) eye(3); zeros(3) zeros(3)];
B = [zeros(3); -inv(J)];
Q = (1/(pi/3)^2)*eye(6); R = eye(size(B,2)); % Q is Tuned Based on Assumed Initial Error
[K,~,~] = lqr(A,B,Q,R);

% Initial Noisy Measurements
rN1 = 2*rand(3,1)-1; rN1 = rN1/norm(rN1);
rN2 = 2*rand(3,1)-1; rN2 = cross(rN1,rN2); rN2 = rN2/norm(rN2);
rN = [rN1,rN2];
[rB_noisy,~] = HSTmeasure([0,0,0,q(:,1)'],V_mag,1,rN); % Magnetometer
[q_noisy,~] = HSTmeasure([0,0,0,q(:,1)'],V_st,2,[]); % Star Tracker
[om_noisy,bias] = HSTmeasure([0,0,0,q(:,1)'],W_gyro,3,[]); % Gyro
    
% MEKF Initial Values
xhist(:,1) = [q(:,1);0;0;0]; % [quaternion; gyro bias]
Phist(:,:,1) = (10*pi/180)^2*eye(6);
whist(:,1) = om_noisy;
yhist(:,1) = [q_noisy; rB_noisy];
W_MEKF = W_gyro;
V_MEKF = blkdiag(V_st, V_mag, V_mag);

% Runge-Kutta Integrator for Controlled Attitude Dynamics
dt = 1; tspan = (0:dt:4*T)/60; % Minutes
x = zeros(length(x0),length(tspan)); x(:,1) = x0; 
theta = zeros(1,length(tspan)); theta(1) = norm(phi0)*180/pi;
u = zeros(3,length(tspan)); xhat = x0(1:6);
for ii = 2:length(tspan)
    % Feedback Control Law
    u(:,ii) = -K*xhat; % u = rho_dot
    % RK4 Step
    k1 = dt*ODEstep([],x(:,ii-1),u(:,ii),A,B,q(:,ii-1),params);
    k2 = dt*ODEstep([],x(:,ii-1)+k1/2,u(:,ii),A,B,q(:,ii-1),params);
    k3 = dt*ODEstep([],x(:,ii-1)+k2/2,u(:,ii),A,B,q(:,ii-1),params);
    k4 = dt*ODEstep([],x(:,ii-1)+k3,u(:,ii),A,B,q(:,ii-1),params);
    x(:,ii) = x(:,ii-1) + (k1+2*k2+2*k3+k4)/6;
    % Compute Angle Error, Quaternion, and Actuator Speeds
    theta(ii) = norm(x(1:3,ii))*180/pi;
    q(:,ii) = qmult(phi2quat(x(1:3,ii)))*qstar;
    inputs(:,ii) = pinv(Bw)*u(:,ii);
    % Noisy Measurements at a Single Time Step
    [rB_noisy,~] = HSTmeasure([0,0,0,q(:,ii)'],V_mag,1,rN); % Magnetometer
    [q_noisy,~] = HSTmeasure([0,0,0,q(:,ii)'],V_st,2,[]); % Star Tracker
    [om_noisy,bias] = HSTmeasure([0,0,0,q(:,ii)'],W_gyro,3,[]); % Gyro

    % MEKF Noisy Measurement Histories
    whist(:,ii) = om_noisy;
    yhist(:,ii) = [q_noisy; rB_noisy];

    % Estimate Attitude with MEKF at a Single Time Step
    [xhist(:,ii), Phist(:,:,ii)] = mekf_step(xhist(:,ii-1), Phist(:,:,ii-1), W_MEKF, V_MEKF, rN, whist(:,ii), yhist(:,ii), dt);
    xhat = [quat2phi(qmult(qconj(qstar))*xhist(1:4,ii)); x(4:6,ii)];
    xhat(1:3,:) = [xhat(2,:); xhat(1,:); -xhat(3,:)];
end
% Check Torque and Speed Values; Rxn Wheel Speed Limit: 3200 RPM(2pi/60) = 335 rad/s
Max_Torque1 = max(inputs(1,:)) % [N-m] Torque Limit per Wheel: 0.82 N-m
Max_Torque2 = max(inputs(2,:))
Max_Torque3 = max(inputs(3,:))
Max_Torque4 = max(inputs(4,:))

% Plot Attitude Error
figure, plot(tspan,theta),ylabel('Degrees'),xlabel('Time (min)')
title('Attitude Error')

% Plot Angular Velocity
om = x(4:6,:); figure
subplot(3,1,1),plot(tspan,om(1,:)),ylabel('\omega_x'),xlabel('Time (s)')
title('HST Angular Velocity with Disturbances')
subplot(3,1,2),plot(tspan,om(2,:)),ylabel('\omega_y'),xlabel('Time (s)')
subplot(3,1,3),plot(tspan,om(3,:)),ylabel('\omega_z'),xlabel('Time (s)')

% Plot Axis-Angle Error
phi = x(1:3,:); figure
subplot(3,1,1),plot(tspan,phi(1,:)),ylabel('\phi_1'),xlabel('Time (s)')
title('Axis-Angle Error')
subplot(3,1,2),plot(tspan,phi(2,:)),ylabel('\phi_2'),xlabel('Time (s)')
subplot(3,1,3),plot(tspan,phi(3,:)),ylabel('\phi_3'),xlabel('Time (s)')
mean((rms(phi,2)))*180/pi

% Plot True vs. Estimated Quaternion
figure,
subplot(4,1,1),plot(tspan,q(1,:),tspan,xhist(1,:),'r-'),xlabel('Time (s)'),ylabel('q_1')
title('HST Quaternion Dynamics'); legend('True', 'Estimated')
subplot(4,1,2),plot(tspan,q(2,:),tspan,xhist(2,:),'r-'),xlabel('Time (s)'),ylabel('q_2')
subplot(4,1,3),plot(tspan,q(3,:),tspan,xhist(3,:),'r-'),xlabel('Time (s)'),ylabel('q_3')
subplot(4,1,4),plot(tspan,q(4,:),tspan,xhist(4,:),'r-'),xlabel('Time (s)'),ylabel('q_4')

% Calculate and Plot Quaternion Estimate Error
e = zeros(3,size(q,2));
for k = 1:size(q,2)
    e(:,k) = quat2phi(qmult(qconj(q(:,k)))*xhist(1:4,k));
end
figure;
subplot(3,1,1);
plot(tspan,(180/pi)*e(1,:)); axis([0 tspan(end) -0.25 0.25]); hold on
plot(tspan,2*(180/pi)*sqrt(squeeze(Phist(1,1,:))),'r');
plot(tspan,-2*(180/pi)*sqrt(squeeze(Phist(1,1,:))),'r');
title('Attitude Error');
subplot(3,1,2);
plot(tspan,(180/pi)*e(2,:)); axis([0 tspan(end) -0.25 0.25]); hold on
plot(tspan,2*(180/pi)*sqrt(squeeze(Phist(2,2,:))),'r');
plot(tspan,-2*(180/pi)*sqrt(squeeze(Phist(2,2,:))),'r');
ylabel('Degrees');
subplot(3,1,3);
plot(tspan,(180/pi)*e(3,:)); axis([0 tspan(end) -0.25 0.25]); hold on
plot(tspan,2*(180/pi)*sqrt(squeeze(Phist(3,3,:))),'r');
plot(tspan,-2*(180/pi)*sqrt(squeeze(Phist(3,3,:))),'r');
xlabel('Time (s)')

% Plot Reaction Wheel Speeds
figure
subplot(4,1,1),plot(tspan,inputs(1,:)),ylabel('\tau_{rxn}^1'),xlabel('t'),title('RXN Wheel Torque')
subplot(4,1,2),plot(tspan,inputs(2,:)),ylabel('\tau_{rxn}^2'),xlabel('t')
subplot(4,1,3),plot(tspan,inputs(3,:)),ylabel('\tau_{rxn}^3'),xlabel('t')
subplot(4,1,4),plot(tspan,inputs(4,:)),ylabel('\tau_{rxn}^4'),xlabel('t')

%% Eigen-Axis Slew
clc,clearvars,close all
load('Parameters\orbital')
load('Parameters\mass')
load('Parameters\geometry')

% Nominal Versine Function
theta_f = 179*pi/180; % Radians
alpha = pi/2000; dt = 0.1; tspan = 0:dt:pi/alpha;
theta = 0.5*theta_f*(1 - cos(alpha*tspan)); % Versine Function
theta_dot = 0.5*alpha*theta_f*sin(alpha*tspan);
theta_ddot = 0.5*(alpha^2)*theta_f*cos(alpha*tspan);

% Nominal Reference Trajectory and Nominal Inputs
r = rand(3,1); r = r/norm(r); % Unit Vector, Axis of Rotation
phi_traj(:,1) = r*theta(1); om_traj(:,1) = r*theta_dot(1);
om_dot_traj(:,1) = r*theta_ddot(1); q_traj(:,1) = phi2quat(phi_traj(:,1));
rho_traj(:,1) = [0; 0; 0]; inputs_traj(:,1) = [0; 0; 0; 0]; speeds_traj(:,1) = [0; 0; 0; 0];
for ii = 2:length(tspan)
    phi_traj(:,ii) = r*theta(ii);
    om_traj(:,ii) = r*theta_dot(ii);
    om_dot_traj(:,ii) = r*theta_ddot(ii);
    q_traj(:,ii) = phi2quat(phi_traj(:,ii));
    rho_dot_traj(:,ii-1) = -J*om_dot_traj(:,ii-1) - cross(om_traj(:,ii-1),J*om_traj(:,ii-1) + rho_traj(:,ii-1));
    rho_traj(:,ii) = rho_traj(:,ii-1) + dt*rho_dot_traj(:,ii-1);

    % Map into Wheel Torque and Wheel Speed
    inputs_traj(:,ii-1) = pinv(Bw)*rho_dot_traj(:,ii-1); % Inverse Dynamics for Wheel Torque
    speeds_traj(:,ii) = pinv(Bw)*rho_traj(:,ii); % Inverse Dynamics for Wheel Speed
end

% Linearize and Discretize About Nominal Trajectory for LQR Gains
Q1 = ((1*pi/180)^-2)*ones(1,3); Q2 = ((1*pi/180)^-2)*ones(1,3); Q3 = ((1*pi/180)^-2)*ones(1,3);
Q = diag([Q1 Q2 Q3]); S = 10*Q;
R = diag([0.8 0.8 0.8].^-2);

for ii = 1:length(tspan)-1 % Integrate "Backwards"
    Ac = [-hat(om_traj(:,ii)) eye(3) zeros(3);
          zeros(3) -inv(J)*(hat(om_traj(:,ii))*J-hat(J*om_traj(:,ii)+rho_traj(:,ii))) -inv(J)*hat(om_traj(:,ii));
          zeros(3) zeros(3) zeros(3)];
    Bc = [zeros(3); -inv(J); eye(3)]; % u = rho_dot
    At = expm(Ac*dt); Bt = Bc*dt;
    [K(:,:,ii),S] = tvlqr(At,Bt,Q,R,S); % Compute Time-History of LQR Gains
end
K = flip(K,3); % Fix Order of Gains
save('nom_traj','phi_traj','om_traj','q_traj','rho_traj','rho_dot_traj','inputs_traj','speeds_traj','K')
disp('Saved Parameters.')

%% Time-Varying LQR for Attitude Tracking with MEKF
clc,clearvars,close all
load('Parameters\orbital')
load('Parameters\mass')
load('Parameters\sensor')
load('Parameters\geometry')
load('nom_traj') % Nominal Values: phi_traj, om_traj, q_traj, rho, rho_dot

% Initial State and Perturbation Conditions
alpha = pi/2000; dt = 0.1; tspan = 0:dt:pi/alpha;
[r_eci, v_eci] = OE2ECI(a, e, i, RAAN, w, anom, mu);
dphi(:,1) = [0; 0; 0];%5*pi/180]; % Normalize and Set Initial Axis-Angle Perturbation
dom(:,1) = [0;0;0]; % rad/s, Zero Initial Angular Velocity Perturbation
drho(:,1) = [0;0;0]; % Zero Initial Wheel Speed Perturbation
dq(:,1) = phi2quat(dphi(:,1)); % Initial Quaternion Perturbation
dx(:,1) = [dphi(:,1); dom(:,1); drho(:,1)]; 

phi(:,1) = phi_traj(:,1) + dphi(:,1);
om(:,1) = om_traj(:,1) + dom(:,1);
rho(:,1) = rho_traj(:,1) + drho(:,1);
q(:,1) = qmult(q_traj(:,1))*dq(:,1);
x(:,1) = [q(:,1); om(:,1); rho(:,1); r_eci; v_eci]; % [q; om; rho; r; v]

% Initial Noisy Measurements
rN1 = 2*rand(3,1)-1; rN1 = rN1/norm(rN1);
rN2 = 2*rand(3,1)-1; rN2 = cross(rN1,rN2); rN2 = rN2/norm(rN2);
rN = [rN1,rN2];
[rB_noisy,~] = HSTmeasure([0,0,0,q(:,1)'],V_mag,1,rN); % Magnetometer
[q_noisy,~] = HSTmeasure([0,0,0,q(:,1)'],V_st,2,[]); % Star Tracker
[om_noisy,bias] = HSTmeasure([0,0,0,q(:,1)'],W_gyro,3,[]); % Gyro
    
% MEKF Initial Values
xhist(:,1) = [q(:,1);0;0;0]; % [quaternion; gyro bias]
Phist(:,:,1) = (10*pi/180)^2*eye(6);
whist(:,1) = om_noisy;
yhist(:,1) = [q_noisy; rB_noisy];
W_MEKF = W_gyro;
V_MEKF = blkdiag(V_st, V_mag, V_mag);

% Runge-Kutta Integrator for Controlled Attitude Tracking
for ii = 2:length(tspan)        
    % Compute Gains and Feedback Control Law
    du = -K(:,:,ii-1)*dx(:,ii-1);
    u(:,ii-1) = rho_dot_traj(:,ii-1) + du; % rho_dot = nominal control 
    
    % RK4 Step
    k1 = ODEtrack([],x(:,ii-1),u(:,ii-1),params);
    k2 = ODEtrack([],x(:,ii-1)+k1/2,u(:,ii-1),params);
    k3 = ODEtrack([],x(:,ii-1)+k2/2,u(:,ii-1),params);
    k4 = ODEtrack([],x(:,ii-1)+k3,u(:,ii-1),params);
    x(:,ii) = x(:,ii-1) + dt*(k1+2*k2+2*k3+k4)/6; % Continuous Orbit Dynamics
    
    % Parse Propagated State [q; om; rho]
    q(:,ii) = x(1:4,ii)/norm(x(1:4,ii));
    phi(:,ii) = quat2phi(q(:,ii));  
    om(:,ii) = x(5:7,ii);
    rho(:,ii) = x(8:10,ii);
    
    % Compute Rotation Angle and Actuator Torques/Speeds
    angle(ii) = norm(phi(:,ii))*180/pi;
    inputs(:,ii-1) = pinv(Bw)*u(:,ii-1); % Scalar Torques
    speeds(:,ii) = pinv(Bw)*x(8:10,ii); % Scalar Speeds
    
    % Noisy Measurements at a Single Time Step
    [rB_noisy,~] = HSTmeasure([0,0,0,q(:,ii)'],V_mag,1,rN); % Magnetometer
    [q_noisy,~] = HSTmeasure([0,0,0,q(:,ii)'],V_st,2,[]); % Star Tracker
    [om_noisy,bias] = HSTmeasure([0,0,0,q(:,ii)'],W_gyro,3,[]); % Gyro
    whist(:,ii) = om_noisy;
    yhist(:,ii) = [q_noisy; rB_noisy];
    
    % Estimate Attitude with MEKF at a Single Time Step
    [xhist(:,ii),Phist(:,:,ii)] = mekf_step(xhist(:,ii-1),Phist(:,:,ii-1),W_MEKF,V_MEKF,rN,whist(:,ii),yhist(:,ii),dt);
    dx(:,ii) = [quat2phi(xhist(1:4,ii)) - phi_traj(:,ii); x(5:7,ii) - om_traj(:,ii); x(8:10,ii) - rho_traj(:,ii)];
end

% Compute Max Torques and Speeds
Max_Torque1 = max(abs(inputs(1,:))) % [N-m] Torque Limit for Each Rxn Wheels: 0.82 N-m
Max_Torque2 = max(abs(inputs(2,:)))
Max_Torque3 = max(abs(inputs(3,:)))
Max_Torque4 = max(abs(inputs(4,:)))

Max_Speed1 = max(abs(speeds(1,:))) % [rad/s] Rxn Wheel Speed Limit: 3200 RPM(2pi/60) = 335 rad/s
Max_Speed2 = max(abs(speeds(2,:)))
Max_Speed3 = max(abs(speeds(3,:)))
Max_Speed4 = max(abs(speeds(4,:)))

% Plot Nominal vs. Closed-Loop Results
% Omega
figure
subplot(3,1,1),plot(tspan,om_traj(1,:),tspan,om(1,:)),xlabel('Time (s)'),ylabel('\omega_x')
title('Nominal vs. Closed-Loop Angular Velocity')
subplot(3,1,2),plot(tspan,om_traj(2,:),tspan,om(2,:)),xlabel('Time (s)'),ylabel('\omega_y')
subplot(3,1,3),plot(tspan,om_traj(3,:),tspan,om(3,:)),xlabel('Time (s)'),ylabel('\omega_z')
legend('Nominal','Closed-Loop','Location','Best')
% Quaternion
figure
subplot(4,1,1),plot(tspan,q_traj(1,:),tspan,q(1,:)),xlabel('Time (s)'),ylabel('q_1')
title('Nominal vs. Closed-Loop Quaternion')
subplot(4,1,2),plot(tspan,q_traj(2,:),tspan,q(2,:)),xlabel('Time (s)'),ylabel('q_2')
subplot(4,1,3),plot(tspan,q_traj(3,:),tspan,q(3,:)),xlabel('Time (s)'),ylabel('q_3')
subplot(4,1,4),plot(tspan,q_traj(4,:),tspan,q(4,:)),xlabel('Time (s)'),ylabel('q_4')
legend('Nominal','Closed-Loop','Location','Best')
% Rotation Angle Theta
figure
plot(tspan,vecnorm(phi_traj)*180/pi,tspan,angle),xlabel('Time (s)'),ylabel('\theta')
title('Nominal vs. Closed-Loop Angle of Rotation')
legend('Nominal','Closed-Loop','Location','Best')
% Axis-Angle
figure
subplot(3,1,1),plot(tspan,phi_traj(1,:),tspan,phi(1,:)),xlabel('Time (s)'),ylabel('\phi_1')
title('Nominal vs. Closed-Loop Axis-Angle')
subplot(3,1,2),plot(tspan,phi_traj(2,:),tspan,phi(2,:)),xlabel('Time (s)'),ylabel('\phi_2')
subplot(3,1,3),plot(tspan,phi_traj(3,:),tspan,phi(3,:)),xlabel('Time (s)'),ylabel('\phi_3')
legend('Nominal','Closed-Loop','Location','Best')
% Wheel Torques
figure, t = tspan(1:end-1);
subplot(4,1,1),plot(t,inputs_traj(1,:),t,inputs(1,:)),xlabel('Time (s)'),ylabel('\tau_1')
title('Nominal vs. Closed-Loop Input Torques')
subplot(4,1,2),plot(t,inputs_traj(2,:),t,inputs(2,:)),xlabel('Time (s)'),ylabel('\tau_2')
subplot(4,1,3),plot(t,inputs_traj(3,:),t,inputs(3,:)),xlabel('Time (s)'),ylabel('\tau_3')
subplot(4,1,4),plot(t,inputs_traj(4,:),t,inputs(4,:)),xlabel('Time (s)'),ylabel('\tau_4')
legend('Nominal','Closed-Loop','Location','Best')