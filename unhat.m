function v = unhat(M)
% Unhat Operator Function
v(1) = -M(2,3);
v(2) = M(1,3);
v(3) = -M(1,2);
end