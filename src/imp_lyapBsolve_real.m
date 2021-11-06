function  U  = imp_lyapBsolve_real( A,B,E,k )

%function X = imp_lyap_solve(A,B);
% 
% Solve  A' X E + E' X A + B'B = 0
%
% Block Hammarling method real variant für implizite Lyapunov Gleichung

n = size(A,1);
m = size(B,1);

[AA,EE,Q,Z]=qz(full(A),full(E),'real');

BB = B*Z;

U = zeros(n,n);

i=1;
k2 = k;

AA_T = AA';
EE_T = EE';

f = zeros(1,n);
qq = zeros(m,n*m);
Z1 = zeros(2,n*2);
Z2 = zeros(m,n*2);
while i<=n
    l = min(n,i+k-1);
    if ~(l==n || (l<n && abs(EE(l+1,l))<10*eps*max(abs(EE(l,l)),abs(EE(l+1,l+1))) && abs(AA(l+1,l))<10*eps*max(abs(AA(l,l)),abs(AA(l+1,l+1)))))
        if (k<k2 || k==1)
            k = k+1;
        else
            k = k-1;
        end
        l = min(n,i+k-1);
    end
    j = i;
    while j<=l
        if (j==l || (j<l && abs(EE(j+1,j))<10*eps*max(abs(EE(j,j)),abs(EE(j+1,j+1))) && abs(AA(j+1,j))<10*eps*max(abs(AA(j,j)),abs(AA(j+1,j+1)))))
            [q,b] = qr(BB(:,j));
            qq(:,(j-1)*m+1:j*m) = q';
            BB(:,j:l) = q'*BB(:,j:l);
            U(j,j) = abs(BB(1,j))/sqrt(-2*(AA(j,j)*EE(j,j)));
            f(j) = BB(1,j)/U(j,j);
            r = -f(j)'*BB(1,j+1:l)-U(j,j)*(EE(j,j)*AA(j,j+1:l)+AA(j,j)*EE(j,j+1:l));
            U(j,j+1:l) = (r/(EE(j,j)*AA(j+1:l,j+1:l)+AA(j,j)*EE(j+1:l,j+1:l)));
            v = U(j,j+1:l)*AA(j+1:l,j+1:l)+AA(j,j+1:l)*U(j,j);
            w = U(j,j+1:l)*EE(j+1:l,j+1:l)+EE(j,j+1:l)*U(j,j);
            BB(1,j+1:l) = (EE(j,j)*v-AA(j,j)*w)/abs(f(j));
            j = j+1;
        else
            [q,b] = qr(BB(:,j:min(j+m-1,n)));
            qq(:,(j-1)*m+1:j*m) = q';
            BB(:,j:l) = q'*BB(:,j:l);
            bb = BB(:,j:j+1)'*BB(:,j:j+1);
            M1 = AA_T(j:(j+1),j:(j+1))*EE(j,j)+EE_T(j:(j+1),j:(j+1))*AA(j,j);
            M2 = AA_T(j:(j+1),j:(j+1))*EE(j+1,j)+EE_T(j:(j+1),j:(j+1))*AA(j+1,j);
            M3 = AA_T(j:(j+1),j:(j+1))*EE(j,j+1)+EE_T(j:(j+1),j:(j+1))*AA(j,j+1);
            M4 = AA_T(j:(j+1),j:(j+1))*EE(j+1,j+1)+EE_T(j:(j+1),j:(j+1))*AA(j+1,j+1);
            M = M2/M4;
            x1 = (M1-M*M3)\(-bb(:,1)+M*bb(:,2));
            x2 = M2\(-bb(:,1)-M1*x1);
            U(j,j) = sqrt(x1(1));
            U(j,j+1) = x1(2)/U(j,j);
            U(j+1,j+1) = sqrt(x2(2)-U(j,j+1)^2);
            Z = U(j:j+1,j:j+1)*EE(j:j+1,j:j+1);
            Z2(:,2*j-1:j*2) = BB(:,j:j+1)/Z;
            Z = U(j:j+1,j:j+1)*AA(j:j+1,j:j+1)/Z;
            Z1(:,2*j-1:j*2) = Z;
            r = -BB(:,j+2:l)'*Z2(:,2*j-1:j*2)-AA_T(j+2:l,j:j+1)*U(j:j+1,j:j+1)'-EE_T(j+2:l,j:j+1)*U(j:j+1,j:j+1)'*Z;
            N1 = AA_T(j+2:l,j+2:l)+EE_T(j+2:l,j+2:l)*Z(1,1);
            N2 = EE_T(j+2:l,j+2:l)*Z(2,1);
            N3 = EE_T(j+2:l,j+2:l)*Z(1,2);
            N4 = AA_T(j+2:l,j+2:l)+EE_T(j+2:l,j+2:l)*Z(2,2);
            N = N2/N4;
            U(j,j+2:l) = (N1-N*N3)\(r(:,1)-N*r(:,2));
            U(j+1,j+2:l) = N2\(r(:,1)-N1*U(j,j+2:l)');
            BB(:,j+2:l) = BB(:,j+2:l)-Z2(:,2*j-1:j*2)*(U(j:j+1,j:j+1)*EE(j:j+1,j+2:l)+U(j:j+1,j+2:l)*EE(j+2:l,j+2:l));
            j = j+2;
        end
    end
    if (l<n)
        p_f = l+1;
        p_l = min(n,l+k);
        if ~(p_l==n || (p_l<n && abs(EE(p_l+1,p_l))<10*eps*max(abs(EE(p_l,p_l)),abs(EE(p_l+1,p_l+1))) && abs(AA(p_l+1,p_l))<10*eps*max(abs(AA(p_l,p_l)),abs(AA(p_l+1,p_l+1)))))
            if (k<k2 || k==1)
                k = k+1;
            else
                k = k-1;
            end
            p_l = min(n,l+k);
        end
        UA = U(1:l,1:l)*AA(1:l,p_f:p_l);
        UE = U(1:l,1:l)*EE(1:l,p_f:p_l);
        j = 1;
        while j<=l
            if (j==l || (j<l && abs(EE(j+1,j))<10*eps*max(abs(EE(j,j)),abs(EE(j+1,j+1))) && abs(AA(j+1,j))<10*eps*max(abs(AA(j,j)),abs(AA(j+1,j+1)))))
                BB(:,p_f:p_l) = qq(:,(j-1)*m+1:j*m)*BB(:,p_f:p_l);
                r = -f(j)*BB(1,p_f:p_l)-(EE(j,j)*UA(j,:)+AA(j,j)*UE(j,:));
                U(j,p_f:p_l) = (r/(EE(j,j)*AA(p_f:p_l,p_f:p_l)+AA(j,j)*EE(p_f:p_l,p_f:p_l)));
                v = U(j,j+1:p_l)*AA(j+1:p_l,p_f:p_l)+AA(j,p_f:p_l)*U(j,j);
                w = U(j,j+1:p_l)*EE(j+1:p_l,p_f:p_l)+EE(j,p_f:p_l)*U(j,j);
                BB(1,p_f:p_l) = (EE(j,j)*v-AA(j,j)*w)/abs(f(j));
                j = j+1;
            else
                BB(:,p_f:p_l) = qq(:,(j-1)*m+1:j*m)*BB(:,p_f:p_l);
                Z = Z1(:,2*j-1:j*2);
                r = -BB(:,p_f:p_l)'*Z2(:,2*j-1:j*2)-UA(j:j+1,:)'-UE(j:j+1,:)'*Z;
                N1 = AA_T(p_f:p_l,p_f:p_l)+EE_T(p_f:p_l,p_f:p_l)*Z(1,1);
                N2 = EE_T(p_f:p_l,p_f:p_l)*Z(2,1);
                N3 = EE_T(p_f:p_l,p_f:p_l)*Z(1,2);
                N4 = AA_T(p_f:p_l,p_f:p_l)+EE_T(p_f:p_l,p_f:p_l)*Z(2,2);
                N = N2/N4;
                U(j,p_f:p_l) = (N1-N*N3)\(r(:,1)-N*r(:,2));
                U(j+1,p_f:p_l) = N2\(r(:,1)-N1*U(j,p_f:p_l)');
                BB(:,p_f:p_l) = BB(:,p_f:p_l)-Z2(:,2*j-1:j*2)*(U(j:j+1,j:j+1)*EE(j:j+1,p_f:p_l)+U(j:j+1,j+2:p_l)*EE(j+2:p_l,p_f:p_l));
                j = j+2;
            end
        end
    end
    i = l+1;
end

U = U*Q;

end

