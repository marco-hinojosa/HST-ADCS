clear;
clc;
close all;

load rBrN.mat
load truth.mat

%% TRIAD Solution

r1 = rN(:,1);
r2 = rN(:,2);
t1 = r1;
t2 = cross(r1,r2)/norm(cross(r1,r2));
TN = [t1,t2,cross(t1,t2)/norm(cross(t1,t2))];

r1 = rB(:,1);
r2 = rB(:,2);
t1 = r1;
t2 = cross(r1,r2)/norm(cross(r1,r2));
TB = [t1,t2,cross(t1,t2)/norm(cross(t1,t2))];

Q_NB = TN*TB'
Q_true

Q_error = Q_true'*Q_NB;
e = unhat(logm(Q_error));
th_e = 180/pi*norm(e)

%% SVD Solution

B = zeros(3);

for i=1:size(rB,2)
    B = B + rB(:,i)*rN(:,i)';
end

[U,S,V] = svd(B');
Q_NB_SVD = U*[1,0,0;0,1,0;0,0,det(U)*det(V)]*V';
Q_NB_SVD*Q_NB_SVD';

Q_error = Q_true'*Q_NB_SVD;
e = unhat(logm(Q_error));
th_e = 180/pi*norm(e)

%% Davenport q-Method Solution

B = zeros(3);
Z = zeros(3,1);
for i=1:size(rB,2)
    B = B + rB(:,i)*rN(:,i)';
    Z = Z + cross(rB(:,i),rN(:,i));
end

K = [B + B' - trace(B)*eye(3) Z;Z' trace(B)]
[eigvec, eigval] = eig(K)
[eigmax, idx] = max(diag(eigval))

q_NB = eigvec(:,idx)/norm(eigvec(:,idx))
v = q_NB(1:3);
s = q_NB(4);

Q_NB_DAVEN = eye(3) + 2*hat(v)*(hat(v) + s*eye(3))
Q_error = Q_true'*Q_NB_DAVEN;
e = unhat(logm(Q_error));
th_e = 180/pi*norm(e)


%% Convex Optimization

cvx_begin sdp
variable Q_cvx(3,3);
maximize trace(B*Q_cvx);
[eye(3), Q_cvx;
    Q_cvx', eye(3)] >= 0;
cvx_end

Q_error = Q_true'*Q_cvx;
e = unhat(logm(Q_error));
th_e = 180/pi*norm(e)

%% Helper Functions

function v = unhat(M)

v(1) = -M(2,3);
v(2) = M(1,3);
v(3) = -M(1,2);

end
function M = hat(v)

M = [0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];

end

