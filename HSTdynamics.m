function x_dot = HSTdynamics(t,init_state,mu,J,rho,params,flag)
% HSTdynamics: Contains full ODE dynamics for all spacecraft states.
%
% Inputs:
%      t - time [s]
%      init_state - initial state vector [r; v; omega; q]
%      mu - central body gravitational parameters [km^3/s^2]
%      J - spacecraft inertia matrix
%      rho - gyro momentum
%
% Outputs:
%   x_dot - linear equation containing ODEs for all states

if length(init_state) == 3 % Angular Velocity
    % Unpack Initial State
    om0 = init_state;
    
    % Matrix Linear Equation
    om_dot = -J\cross(om0,J*om0 + rho);
    x_dot = om_dot;
elseif length(init_state) == 6 % Linear Position and Velocity
    % Unpack Initial State
    rvec = init_state(1:3);
    vvec = init_state(4:6);   
    r = norm(rvec);
    
    % Matrix Linear Equation
    rv_dot = [zeros(3,3) eye(3,3);
        -(mu/(r^3))*eye(3,3) zeros(3,3)]*[rvec; vvec];
    x_dot = rv_dot;  
elseif length(init_state) == 7 % Angular Velocity and Quaternion
    % Unpack Initial State
    om0 = init_state(1:3);
    q0 = init_state(4:7);
    qhat = qmult(q0);
    
    % Matrix Linear Equation
    om_dot = -J\cross(om0,J*om0 + rho);
    q_dot = (1/2)*qhat*[om0; 0];
    x_dot = [om_dot; q_dot];
elseif length(init_state) == 13 % All States
    % Unpack Initial State
    om0 = init_state(1:3);
    q0 = init_state(4:7);
    rvec = init_state(8:10);
    vvec = init_state(11:13);   
    qhat = qmult(q0);
    r = norm(rvec);
    
    if flag == 1
        % Compute Gravity Gradient Torque
        Q = quat2Q(q0); % nQb
        Tau_g = cross(3*mu/(dot(rvec,rvec)^(5/2))*rvec, Q*J/(1000^2)*rvec); % [N-km]
        Tau_g = Q'*Tau_g*1000; % Rotate into Body Frame and Convert Units, [N-m]
        
        S = params(1,:);
        n = params(2:4,:);
        c = params(5:7,:);
        [D,Tau_D] = HSTdrag(rvec,vvec,Q,S,n,c);
    else
        Tau_g = zeros(3,1);
        D = zeros(3,1);
        Tau_D = zeros(3,1);
    end
    D = Q*D/1000; % Rotate Drag into Inertial Frame, [kg-km/s^2]
    M = 11100; % HST Mass, [kg]
    
    % Matrix Linear Equation
    om_dot = J\(Tau_g + Tau_D - cross(om0,J*om0 + rho));
    q_dot = (1/2)*qhat*[om0; 0];
    rv_dot = [zeros(3,3) eye(3,3);
        -(mu/(r^3))*eye(3,3) zeros(3,3)]*[rvec; vvec] - [zeros(3,3); eye(3,3)]*(D/M);
    x_dot = [om_dot; q_dot; rv_dot];
end
t
end
