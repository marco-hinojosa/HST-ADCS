function [y,C] = measurement(q,rN)

Q = quat2Q(q);
rB = [];
for ii = 1:size(rN,2)
    rB = [rB; Q'*rN(:,ii)];
end
C = [eye(3), zeros(3)];
for ii = 1:size(rN,2)
    C = [C; hat(rB(3*ii-2:3*ii)) zeros(3)];
end
y = rB(:);
end