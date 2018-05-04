function x_dot = HSTdynamics(mu,J,rho,init_state)

% Define State x = [r; r_dot; om; quat], x_dot = [r_dot; r_ddot; om_dot; quat_dot]
rvec = init_state(1:3);
vvec = init_state(4:6);
om0 = init_state(7:9);
q0 = init_state(10:13);

v = q0(1:3);
s = q0(4);
qhat = [s*eye(3)+hat(v) v;-v' s];
r = norm(rvec);

% Matrix Linear Equation
x_dot = [zeros(3,3) eye(3,3);
    -(mu/(r^3))*eye(3,3) zeros(3,3)]*[rvec; vvec];

x_dot(7:9,1) = -J\cross(om0,J*om0 + rho);
x_dot(10:13,1) = (1/2)*qhat*[om0; 0];