function x_dot = HSTorbit(mu,init_state)

% Define State x = [r; r_dot], x_dot = [r_dot; r_ddot]
rvec = init_state(1:3);
vvec = init_state(4:6);
r = norm(rvec);

% Matrix Linear Equation
x_dot = [zeros(3,3) eye(3,3);
    -(mu/(r^3))*eye(3,3) zeros(3,3)]*[rvec; vvec];

% State Equations
%x1 = rvec;
%x2 = vvec;

%x1_dot = x2; % r_dot (= vvec)
%x2_dot = -mu*x1/(r^3); % r_ddot

%x_dot = [x1_dot; x2_dot];