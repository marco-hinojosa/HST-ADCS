function om_dot = HSTspin(t,om,J,rho)
% Matrix Linear Equation
om_dot = -J\cross(om,J*om + rho);
