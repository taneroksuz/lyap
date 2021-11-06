function [ X,res ] = imp_lyap1solve( A,B,E )

%function X = imp_lyap_solve(A,B);
% 
% Solve  A' X E + E' X A + B = 0
%
% 1-Solve method fuer obige implizite Lyapunov Gleichung

n = size(A,1); 

[AA,EE,Q,Z]=qz(A,E,'complex');

BB = Z'*B*Z;

AA_T = AA';
EE_T = EE';

X = zeros(n,n);

for k=1:n
    d = AA_T*X(:,1:k-1)*EE(1:k-1,k)+EE_T*X(:,1:k-1)*AA(1:k-1,k);
    X(:,k)=(EE(k,k)*AA_T+AA(k,k)*EE_T)\(-BB(:,k)-d);
end

X = Q'*X*Q;
res = norm(A'*X*E+E'*X*A+B,'fro');


end

