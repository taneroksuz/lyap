function [X,res] = imp_lyap_solve( A,B,E )

%function X = imp_lyap_solve(A,B);
% 
% Solve  A' X E + E' X A + B = 0
%
% Bartels Stewart method für implizite Lyapunov Gleichung

n = size(A,1); 

[AA,EE,Q,Z]=qz(A,E,'complex');

BB = Z'*B*Z;

X = zeros(n,n);

for i=1:n
    b_11 = BB(i,i);
    e_11 = EE(i,i);
    a_11 = AA(i,i);
    b = BB(i,i+1:n)';
    e = EE(i,i+1:n)';
    a = AA(i,i+1:n)';
    B_1 = BB(i+1:n,i+1:n);
    E_1 = EE(i+1:n,i+1:n);
    A_1 = AA(i+1:n,i+1:n);
    x_11 = -b_11/(a_11*e_11'+e_11*a_11');
    x = (e_11*A_1'+a_11*E_1')\(-b-x_11*(a_11*e+e_11*a));
    B_1 = B_1+x_11*(a*e'+e*a')+(A_1'*x*e'+e*x'*A_1)+(E_1'*x*a'+a*x'*E_1);
    BB(i+1:n,i+1:n) = B_1;
    X(i,i) = x_11;
    X(i,i+1:n) = x';
    X(i+1:n,i) = x;
end

X = Q'*X*Q;
res = norm(A'*X*E+E'*X*A+B,'fro');
end

