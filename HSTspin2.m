function om_dot = HSTspin2(t,om,J,rho)
% Matrix Linear Equation
om_dot = -J\cross(om,J*om + rho);

