function [ X,res ] = imp_lyap2solve_real( A,B,E )

%function X = imp_lyap_solve(A,B);
% 
% Solve  A' X E + E' X A + B = 0
%
% 2-Solve method fuer obige implizite Lyapunov Gleichung

n = size(A,1); 

[AA,EE,Q,Z]=qz(A,E,'real');

BB = Z'*B*Z;

AA_T = AA';
EE_T = EE';

X = zeros(n,n);
k = 1;
while (k<(n+1))
    if (k==n || (k<n && abs(EE(k+1,k))<10*eps*max(abs(EE(k,k)),abs(EE(k+1,k+1))) && abs(AA(k+1,k))<10*eps*max(abs(AA(k,k)),abs(AA(k+1,k+1)))))
        d = AA_T*X(:,1:k-1)*EE(1:k-1,k)+EE_T*X(:,1:k-1)*AA(1:k-1,k);
        X(:,k)=(EE(k,k)*AA_T+AA(k,k)*EE_T)\(-BB(:,k)-d);
        k = k+1;
    else
        d = AA_T*X(:,1:k-1)*EE(1:k-1,k:k+1)+EE_T*X(:,1:k-1)*AA(1:k-1,k:k+1);
        b = -BB(:,k:k+1)-d;
        M1 = AA_T*EE(k,k)+EE_T*AA(k,k);
        M2 = EE_T*AA(k+1,k);
        M3 = EE_T*AA(k,k+1);
        M4 = AA_T*EE(k+1,k+1)+EE_T*AA(k+1,k+1);
        M = M2/M4;
        X(:,k) = (M1-M*M3)\(b(:,1)-M*b(:,2));
        X(:,k+1) = M2\(b(:,1)-M1*X(:,k));
        k = k+2;
    end
end

X = Q'*X*Q;
res = norm(A'*X*E+E'*X*A+B,'fro');


end

