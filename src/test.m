clear all;

n = 400;
m = 2;
k = round(n/5);

A = rand(n);
B = eye(m,n);
e = max(real(eig(A)));
if (e>0)
    A = A-(e+0.001)*eye(n);
end

BB = B'*B;

profile on;
tic;
[U1,Q1] = lyapBsolve(A',B,k);
t1 = toc;


tic;
X2 = lyap2solve(A,BB);
t2 = toc;

r2 = norm(A*X2+X2*A'+BB);


% tic;
% [X3,r3] = lyap_solve(A,BB);
% t3 = toc;

tic;
X4 = lyap(A,BB);
t4 = toc;
r4 = norm(A*X4+X4*A'+BB,'fro');

profile off;

X1 = (U1'*U1);

diff_X = abs(X1-Q1'*X2*Q1);

X1 = Q1*X1*Q1';
r1 = norm(A*X1+X1*A'+B'*B);

profile viewer;
