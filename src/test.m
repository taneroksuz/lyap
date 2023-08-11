clc;
close all;
clear all;
profile off;
profile clear;

pkg load control;
warning('error', 'all');

n = 2000;
m = 1;
k = 100;

A = rand(n);
B = eye(m,n);
e = max(real(eig(A)));
if (e>0)
    A = A-(e+0.001)*eye(n);
end

BB = B'*B;

profile on;

tic;
X1 = lyapBHsolve(A',B,k);
t1 = toc;
r1 = norm(A*X1+X1*A'+BB,'fro');

% tic;
% X2 = lyap2solve(A,BB);
% t2 = toc;
% r2 = norm(A*X2+X2*A'+BB,'fro');

% tic;
% X3 = lyapBSsolve(A,BB);
% t3 = toc;
% r3 = norm(A*X3+X3*A'+BB,'fro');

% tic;
% X4 = lyap(A,BB);
% t4 = toc;
% r4 = norm(A*X4+X4*A'+BB,'fro');

profile off;

profshow();
