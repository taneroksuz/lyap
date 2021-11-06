function X = gen_lyap2solve(A,B,E)

% function X = gen_lyap2solve(A,B,E);
% 
% Solve  A' X E + E' X A + B = 0
%
% 2-solve method for generalized lyapunov equation

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
        [aa,ee,q,z] = qz(AA(k:k+1,k:k+1),EE(k:k+1,k:k+1));
        bb = b*z;
        M1 = AA_T*ee(1,1)+EE_T*aa(1,1);
        M2 = AA_T*ee(1,2)+EE_T*aa(1,2);
        M3 = AA_T*ee(2,2)+EE_T*aa(2,2);
        X(:,k) = M1\bb(:,1);
        X(:,k+1) = M3\(bb(:,2)-M2*X(:,k));
        X(:,k:k+1) = X(:,k:k+1)*q;
        k = k+2;
    end
end

X = Q'*X*Q;

end

