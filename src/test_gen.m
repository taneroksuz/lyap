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

A = randn(n);
E = randn(n);
e = max(real(eig(A,E)));
if (e>0)
    A = A-(2*e)*E;
end

B = eye(m,n);
BB = B'*B;

profile on;

% tic;
% X_lyap1 = gen_lyap1solve(A,BB,E);
% r_lyap1 = norm(A'*X_lyap1*E+E'*X_lyap1*A+BB,'fro');
% t_lyap1 = toc;

% tic;
% X_lyap2 = gen_lyap2solve(A,BB,E);
% r_lyap2 = norm(A'*X_lyap2*E+E'*X_lyap2*A+BB,'fro');
% t_lyap2 = toc;

% tic;
% X_lyap2_real = gen_lyap2solve_real(A,BB,E);
% r_lyap2_real = norm(A'*X_lyap2_real*E+E'*X_lyap2_real*A+BB,'fro');
% t_lyap2_real = toc;

% tic;
% X_lyapBH = gen_lyapBHsolve(A,B,E,k);
% r_lyapBH = norm(A'*X_lyapBH*E+E'*X_lyapBH*A+BB,'fro');
% t_lyapBH = toc;

tic;
X_lyapBH_real = gen_lyapBHsolve_real(A,B,E,k);
r_lyapBH_real = norm(A'*X_lyapBH_real*E+E'*X_lyapBH_real*A+BB,'fro');
t_lyapBH_real = toc;

% tic;
% As = A/E;
% BBs = (E'\(BB/E));
% X_lyap_inv = lyap(As',BBs);
% r_lyap_inv = norm(A'*X_lyap_inv*E+E'*X_lyap_inv*A+BB,'fro');
% t_lyap_inv = toc;

% tic;
% X_lyap = lyap(A',BB',[],E');
% r_lyap = norm(A'*X_lyap*E+E'*X_lyap*A+BB,'fro');
% t_lyap = toc;

profile off;

profshow();
