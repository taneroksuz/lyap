function X = lyapBHsolve(A,B,k)

% function X = lyapBHsolve(A,B,k);
% 
% Solve  A' X + X A + B'B = 0
%
% Kressner Block Hammarling method for lyapunov equation

n = size(A,1);
m = size(B,1);

[Q,AA]=schur(full(A),'complex');

BB = B*Q;

U = zeros(n,n);

f = zeros(1,n);
qq = zeros(m,n*m);
i=1;
while i<=n
    l = min(n,i+k-1);
    for j=i:l
        [q,b] = qr(BB(:,j));
        qq(:,(j-1)*m+1:j*m) = q';
        BB(:,j:l) = q'*BB(:,j:l);
        U(j,j) = BB(1,j)/sqrt(-(AA(j,j)'+AA(j,j)));
        if (abs(U(j,j)) == 0)
            U(j,j) = 0;
            break;
        end
        f(j) = BB(1,j)/U(j,j);
        r = -f(j)'*BB(1,j+1:l)-U(j,j)'*AA(j,j+1:l);
        U(j,j+1:l) = (r/(AA(j+1:l,j+1:l)+AA(j,j)'*eye(l-j)));
        BB(1,j+1:l) = BB(1,j+1:l)-f(j)*U(j,j+1:l);
    end
    if (l<n)
        p_f = l+1;
        p_l = min(n,l+k);
        U(1:l,p_f:p_l) = U(1:l,1:l)*AA(1:l,p_f:p_l);
        for j=1:l
            BB(:,p_f:p_l) = qq(:,(j-1)*m+1:j*m)*BB(:,p_f:p_l);
            rp = -f(j)'*BB(1,p_f:p_l)-U(j,p_f:p_l);
            U(j,p_f:p_l) = (rp/(AA(p_f:p_l,p_f:p_l)+AA(j,j)'*eye(p_l-p_f+1)));
            if (norm(U(j,p_f:p_l)) == 0)
                break;
            end
            BB(1,p_f:p_l) = BB(1,p_f:p_l)-f(j)*U(j,p_f:p_l);
        end
    end
    i = l+1;
end

% U(isnan(U)) = 0;

U = U*Q';
X = U'*U;

end

