clc;
close all;
clear all;

n = 200;
% m = 2;

A = randn(n);
E = randn(n);
e = max(real(eig(A,E)));
if (e>0)
    A = A-(2*e)*E;
end


X = randn(n);
X = X'+X;
x = min(eig(X));
if (x<0)
    X = X-(2*x)*eye(n);
end

% B = eye(m,n);
% BB = B'*B;
BB = -A'*X*E;
BB = BB+BB';

profile on;

tic;
[X_lyap1,r_lyap1] = imp_lyap1solve(A,BB,E);
e_lyap1 = norm(X_lyap1-X,'fro');
t_lyap1 = toc;

tic;
[X_lyap2,r_lyap2] = imp_lyap2solve(A,BB,E);
e_lyap2 = norm(X_lyap2-X,'fro');
t_lyap2 = toc;

tic;
[X_lyap2_real,r_lyap2_real] = imp_lyap2solve_real(A,BB,E);
e_lyap2_real = norm(X_lyap2_real-X,'fro');
t_lyap2_real = toc;

tic;
[X_lyap2_real_2,r_lyap2_real_2] = imp_lyap2solve_real_2(A,BB,E);
e_lyap2_real_2 = norm(X_lyap2_real_2-X,'fro');
t_lyap2_real_2 = toc;

tic;
As = A/E;
BBs = (E'\(BB/E));
X_lyap_inv = lyap(As',BBs);
r_lyap_inv = norm(A'*X_lyap_inv*E+E'*X_lyap_inv*A+BB,'fro');
e_lyap_inv = norm(X_lyap_inv-X,'fro');
t_lyap_inv = toc;

tic;
X_lyap = lyap(A',BB',[],E');
r_lyap = norm(A'*X_lyap*E+E'*X_lyap*A+BB,'fro');
e_lyap = norm(X_lyap-X,'fro');
t_lyap = toc;

profile off;
profile viewer;
