%% Reference Eigen-Axis Slew and Controls for 180° Maneuver
clc, clearvars, close all
load mass.mat
% Declare Slew Angle and Initial Conditions
qd = [0 0 0 1]'; % Desired Attitude
r = rand(3,1); r = r/norm(r); % Random Axis of Rotation
theta0 = pi*.99; % 180° Maneuver
phi0 = theta0*r; 
q0 = phi2quat(phi0); q0 = q0/norm(q0);

% Versine Function
tf = 1800; dtOpt = .2; t = 0:dtOpt:tf; alpha = pi/tf; % Time Parameters
thetaOpt = theta0 - theta0*1/2*(1-cos(alpha*t));
phiOpt = thetaOpt.*r;

% Create Eigen-Axis Slew
omOpt = zeros(3,length(t)); qOpt = zeros(4,length(t)); rhoOpt = zeros(3,length(t));
for ii = 1:length(t)
    qOpt(:,ii) = phi2quat(phiOpt(:,ii));
    omOpt(:,ii) = -0.5*r*theta0*alpha*sin(alpha*t(ii));
    domOpt(:,ii) = -0.5*r*theta0*alpha*alpha*cos(alpha*t(ii));
    tauOpt(:,ii) = -J*domOpt(:,ii) - cross(omOpt(:,ii), J*omOpt(:,ii) + rhoOpt(:,ii));
    rhoOpt(:,ii+1) = rhoOpt(:,ii) + dtOpt*tauOpt(:,ii);
end
rhoOpt = rhoOpt(:,1:end-1); % Chop off the extra final value
save('TrajOptParams','J','dtOpt','q0','qd','qOpt','omOpt','rhoOpt','tauOpt');
disp('Parameters Saved.')

%% Attitude Trajectory Optimization with Reaction Wheels
close all
load('TrajOptParams.mat')

% Discretize at Interpolated Points
q = qOpt(:,1:90:end);
om = omOpt(:,1:90:end);
rho = rhoOpt(:,1:90:end);
tau = tauOpt(:,1:90:end);
dt = dtOpt*90;

% Initialize Number of Collocation Points
N = size(q,2);
Nx = N*14 - 4;
Nc = 10*(N-1)+N-2;

% Initial Guess
x0 = [q; om; rho; tau; dt*ones(1,N)]; x0 = x0(:);
x0(end-3:end) = []; 

% Upper and Lower Bounds on the State: [q; om; rho; tau; dt] 14x1
lb = [-Inf*ones(4,N); -100*ones(3,N); -4000*ones(3,N); -0.8*ones(3,N); dt/10*ones(1,N)];
ub = [Inf*ones(4,N); 100*ones(3,N); 4000*ones(3,N); 0.8*ones(3,N); dt*ones(1,N)];
lb = lb(:); ub = ub(:);
lb(end-3:end) = []; ub(end-3:end) = [];

% State Equality Constraints at Boundaries and Knot Points
Aeq = zeros(N-2,14*N-4);
for ii = 1:N-2 % Constraint on Time Steps
    Aeq(ii,14*ii) = 1;
    Aeq(ii,14*(ii+1)) = -1;
end
beq = zeros(N-2,1); % Constraint on Time Step
Aeq = [eye(10), zeros(10,14*(N-1)); zeros(10,14*(N-1)) eye(10); Aeq]; % Initial and Final BC
beq = [x0(1:10); q(:,end); om(:,end); rho(:,end); beq];  % Initial and Final BC

cl = zeros(Nc,1);
cu = zeros(Nc,1);
for ii = 11:11:Nc
    cl(ii) = 1;
    cu(ii) = 1;
end
    
objGrad = zeros(1,14*N-4);
for ii = 1:length(objGrad)
    if mod(ii,14) == 0
        objGrad(ii) = 1;
    end
end

% Build and Solve OPTI Problem
opts = optiset('solver','ipopt','display','iter');
Opt = opti('fun',@(x) fun(x,objGrad),'grad',@(x) grad(x,objGrad),'nl',@(x) nl(x,N,J,q0,qd), ...
               cl,cu,'jac',@(x) jac(x,N,J,q0,qd),'bounds',lb,ub,'x0',x0, ...
               'eq',Aeq,beq,'options',opts);
        
[x,fval,exitflag,info] = solve(Opt);     

% Parse the Optimal Solution and Save Results
xmod = reshape(x(1:14*(N-1)),14,N-1);
qs = [xmod(1:4,:) x(end-9:end-6)];
ws = [xmod(5:7,:) x(end-5:end-3)];
rhos = [xmod(8:10,:) x(end-2:end)];
us = xmod(11:13,:);
dts = xmod(14,:);
ts = 0:dts(1):sum(dts)+.001;
tOpt = 0:dtOpt:(length(qOpt)-1)*dtOpt;

save('TrajOptResults','qs','ws','rhos','us','ts','qOpt','omOpt','tauOpt','tOpt')

%% Plot Results
load('TrajOptParams.mat')
load('TrajOptResults.mat')
load('geometry.mat')
for ii = 1:length(rhoOpt)
    hOpt(:,ii) = pinv(Bw)*rhoOpt(:,ii);
end
for jj = 1:length(rhos)
    hs(:,jj) = pinv(Bw)*rhos(:,jj);
end

figure,subplot(4,1,1)
plot(tOpt,qOpt(1,:),'LineWidth',2),hold on
plot(ts,qs(1,:),'LineWidth',2),title('Quaternion')
legend('Eigen-Axis Slew','Time-Optimal Slew','Location','Best')
ylabel('q_1')
subplot(4,1,2)
plot(tOpt,qOpt(2,:),'LineWidth',2),hold on
plot(ts,qs(2,:),'LineWidth',2), ylabel('q_2')
subplot(4,1,3)
plot(tOpt,qOpt(3,:),'LineWidth',2),hold on
plot(ts,qs(3,:),'LineWidth',2), ylabel('q_3')
subplot(4,1,4)
plot(tOpt,qOpt(4,:),'LineWidth',2),hold on
plot(ts,qs(4,:),'LineWidth',2), ylabel('q_4')
xlabel('Time (s)')

figure,subplot(3,1,1)
plot(tOpt,omOpt(1,:),'LineWidth',2),hold on
plot(ts,ws(1,:),'LineWidth',2),title('Omega')
legend('Eigen-Axis Slew','Time-Optimal Slew','Location','Best')
subplot(3,1,2)
plot(tOpt,omOpt(2,:),'LineWidth',2),hold on
plot(ts,ws(2,:),'LineWidth',2)
ylabel('Angular Velocity (rad/s)')
subplot(3,1,3)
plot(tOpt,omOpt(3,:),'LineWidth',2),hold on
plot(ts,ws(3,:),'LineWidth',2)
xlabel('Time (s)')

figure,subplot(3,1,1)
plot(tOpt,tauOpt(1,:),'LineWidth',2),hold on
plot(ts(1:end-1),us(1,:),'LineWidth',2),title('Control')
legend('Eigen-Axis Slew','Time-Optimal Slew','Location','Best')
subplot(3,1,2)
plot(tOpt,tauOpt(2,:),'LineWidth',2),hold on
plot(ts(1:end-1),us(2,:),'LineWidth',2)
ylabel('Torque (N-m)')
subplot(3,1,3)
plot(tOpt,tauOpt(3,:),'LineWidth',2),hold on
plot(ts(1:end-1),us(3,:),'LineWidth',2)
xlabel('Time (s)')

figure,subplot(4,1,1)
plot(tOpt,hOpt(1,:),'LineWidth',2),hold on
plot(ts,hs(1,:),'LineWidth',2),title('Reaction Wheel Momentum')
legend('Eigen-Axis Slew','Time-Optimal Slew','Location','Best'), ylabel('h_1')
subplot(4,1,2)
plot(tOpt,hOpt(2,:),'LineWidth',2),hold on
plot(ts,hs(2,:),'LineWidth',2), ylabel('h_2')
subplot(4,1,3)
plot(tOpt,hOpt(3,:),'LineWidth',2),hold on
plot(ts,hs(3,:),'LineWidth',2), ylabel('h_3')
subplot(4,1,4)
plot(tOpt,hOpt(4,:),'LineWidth',2),hold on
plot(ts,hs(4,:),'LineWidth',2), ylabel('h_4')
xlabel('Time (s)')

%% Auxiliary Functions
function obj = fun(x,objGrad)
obj = objGrad*x;
end

function objGrad = grad(x,objGrad)
objGrad = objGrad;
end

% Compute Nonlinear Constraints
function [ceq] = nl(x,N,J,q0,qd)
[ceq, ~] = nonlcon(x,N,J,q0,qd);
end

% Compute Nonlinear Constraint Jacobians
function [dceq] = jac(x,N,J,q0,qd)
[~, dceq] = nonlcon(x,N,J,q0,qd);
end

% Define Nonlinear Constraints
function [ceq, dceq] = nonlcon(x,N,J,q0,qd)
Nc = 10*(N-1)+N-2;
xmod = reshape(x(1:14*(N-1)),14,N-1);
q = [xmod(1:4,:) x(end-9:end-6)];
w = [xmod(5:7,:) x(end-5:end-3)];
rho = [xmod(8:10,:) x(end-2:end)];
u = xmod(11:13,:);
dt = xmod(14,:);

ceq = zeros(Nc,1);
dceq = zeros(Nc,14*N-4);
for ii = 1:N-1
    qm = .5*(q(:,ii) + q(:,ii+1));
    wm = .5*(w(:,ii) + w(:,ii+1));
    rhom = .5*(rho(:,ii) + rho(:,ii+1));
    [fx, dfdx, dfdu] = dynamics(qm,wm,rhom,u(:,ii),J);
    ceq(1+11*(ii-1):10+11*(ii-1),1) = [q(:,ii);w(:,ii);rho(:,ii)]+fx*dt(ii) - [q(:,ii+1);w(:,ii+1);rho(:,ii+1)];
    if ii ~= (N-1)
        ceq(11*ii,1) = q(:,ii)'*q(:,ii); % Quaternion Norm Constraint
    end
    A = []; A = fillrows(fx, dfdx, dfdu, dt(ii), q(:,ii), ii, length(x));
    if ii ~= (N-1)
        dceq(1+11*(ii-1):11+11*(ii-1),:) =  A;
    else
        dceq(1+11*(ii-1):10+11*(ii-1),:) =  A;
    end
end
end

% System Dynamics
function [fx,dfx,du] = dynamics(q,w,rho,u,J)
fx(1:4,1) = .5*qmult(q)*[w;0];
fx(5:7,1) = -inv(J)*(u + cross(w,J*w + rho));
fx(8:10,1) = u;

dfx = [.5*qmultR([w;0]), .5*qmult(q)*[eye(3); zeros(1,3)],    zeros(4,3);
       zeros(3,4),       -inv(J)*(hat(w)*J - hat(J*w + rho)), -inv(J)*hat(w);
       zeros(3,4),       zeros(3),                            zeros(3)];

du = [zeros(4,3); -inv(J); eye(3)];
end

% Create Sparse Matrix
function A = fillrows(fx, dfx, du, dt, q, n, N)
m = n-1;
if n ~= (N+4)/14 - 1
    A = zeros(11,N);
    A(1:10,m*14+1:m*14+10) = eye(10) + dfx*.5*eye(10)*dt;
    A(1:10,m*14+11:m*14+13) = dt*du;
    A(1:10,m*14+14) = fx;
    A(1:10,m*14+15:m*14+24) = dfx*.5*eye(10)*dt - eye(10);
    A(11,m*14+1:m*14+4) = 2*q';  % Quaternion Norm Constraint
else
    A = zeros(10,N);
    A(1:10,m*14+1:m*14+10) = eye(10) + dfx*.5*eye(10)*dt;
    A(1:10,m*14+11:m*14+13) = dt*du;
    A(1:10,m*14+14) = fx;
    A(1:10,m*14+15:m*14+24) = dfx*.5*eye(10)*dt - eye(10);
end
end
