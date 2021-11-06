function  U  = imp_lyapBsolve( A,B,E,k )

%function X = imp_lyap_solve(A,B);
% 
% Solve  A' X E + E' X A + B'B = 0
%
% Block Hammarling method fuer implizite Lyapunov Gleichung

n = size(A,1);
m = min(size(B));

[AA,EE,Q,Z]=qz(full(A),full(E),'complex');

BB = B*Z;

U = zeros(n,n);

i=1;

f = zeros(1,n);
qq = zeros(m,n*m);
while i<=n
    l = min(n,i+k-1);
    for j=i:l
        [q,b] = qr(BB(:,j));
        qq(:,(j-1)*m+1:j*m) = q';
        BB(:,j:l) = q'*BB(:,j:l);
        U(j,j) = abs(BB(1,j))/sqrt(-(AA(j,j)*EE(j,j)'+AA(j,j)'*EE(j,j)));
        f(j) = BB(1,j)/U(j,j);
        r = -f(j)'*BB(1,j+1:l)-U(j,j)*(EE(j,j)'*AA(j,j+1:l)+AA(j,j)'*EE(j,j+1:l));
        U(j,j+1:l) = (r/(EE(j,j)'*AA(j+1:l,j+1:l)+AA(j,j)'*EE(j+1:l,j+1:l)));
        v = U(j,j+1:l)*AA(j+1:l,j+1:l)+AA(j,j+1:l)*U(j,j)';
        w = U(j,j+1:l)*EE(j+1:l,j+1:l)+EE(j,j+1:l)*U(j,j)';
        BB(1,j+1:l) = (EE(j,j)*v-AA(j,j)*w)/abs(f(j));
    end
    if (l<n)
        p_f = l+1;
        p_l = min(n,l+k);
        UA = U(1:l,1:l)*AA(1:l,p_f:p_l);
        UE = U(1:l,1:l)*EE(1:l,p_f:p_l);
        for j=1:l
            BB(:,p_f:p_l) = qq(:,(j-1)*m+1:j*m)*BB(:,p_f:p_l);
            r = -f(j)'*BB(1,p_f:p_l)-(EE(j,j)'*UA(j,:)+AA(j,j)'*UE(j,:));
            U(j,p_f:p_l) = (r/(EE(j,j)'*AA(p_f:p_l,p_f:p_l)+AA(j,j)'*EE(p_f:p_l,p_f:p_l)));
            v = U(j,j+1:p_l)*AA(j+1:p_l,p_f:p_l)+AA(j,p_f:p_l)*U(j,j)';
            w = U(j,j+1:p_l)*EE(j+1:p_l,p_f:p_l)+EE(j,p_f:p_l)*U(j,j)';
            BB(1,p_f:p_l) = (EE(j,j)*v-AA(j,j)*w)/abs(f(j));
        end
    end
    i = l+1;
end

U = U*Q;

end

