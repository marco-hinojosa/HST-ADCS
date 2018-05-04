function om_dot = HSTspin(t,om,J,rho)
Jvec = diag(J);
% Matrix Linear Equation
om_dot1 = -(Jvec(3) - Jvec(2))*om(2)*om(3)/Jvec(1);
om_dot2 = -(Jvec(1) - Jvec(3))*om(1)*om(3)/Jvec(2);
om_dot3 = -(Jvec(2) - Jvec(1))*om(1)*om(2)/Jvec(3);

om_dot = [om_dot1;om_dot2;om_dot3];

