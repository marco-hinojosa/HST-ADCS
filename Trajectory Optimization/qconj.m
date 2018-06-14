function q = qconj(qq)
% Quaternion Conjugate Operator
v = qq(1:3);
v = -v;
q = [v; qq(4)];
end
