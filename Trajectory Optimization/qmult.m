function qhat = qmult(q)
% Quaternion Hat Operator
% Convert Quaternion into a 4x4 Square Matrix
v = q(1:3);
s = q(4);
qhat = [s*eye(3) + hat(v)  v ; -v' s]; 
end