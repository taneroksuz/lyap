function X = lyap2solve(A,B);
%function X = lyap2solve(A,B);
% 
% Solve  A X + X A' + B = 0
%
% Zhou and Sorensen 2-solve  method

%Ech = chol(E,'lower');%Ech*Ech' = E;
%A = (Ech\A)/(Ech'); B = (Ech\B)/(Ech');


m = size(A,1); n = size(B,2);
[Q,R]=schur(A);
idx = [m:-1:1];Q2=Q(:,idx); R2=R(idx,idx)';
B = Q'*B*Q2;

Rsq = R*R; I=speye(m); 
j=1;

while (j < n+1)
 
  if j==n | (j<n  & abs(R2(j+1,j))<10*eps*max( abs(R2(j,j)), abs(R2(j+1,j+1))) )
 
      if (j>1), b = -B(:,j) - X(:,1:j-1)*R2(1:j-1,j);else
              b = -B(:,j) ;end
      X(:,j) = (R+R2(j,j)*I)\b;
      j = j +1;
  else

      r11 = R2(j,j);  r12 = R2(j,j+1);
      r21 = R2(j+1,j); r22 = R2(j+1,j+1);

      if (j>1), b = -B(:,j:j+1) - X(:,1:j-1)*R2(1:j-1,j:j+1); else
                b = -B(:,j:j+1); end
      b = [R*b(:,1)+r22*b(:,1)-r21*b(:,2), R*b(:,2)+r11*b(:,2)-r12*b(:,1)];
      X(:,j:j+1) = ( Rsq+(r11+r22)*R + (r11*r22-r12*r21)*I)\b;
      j = j + 2;
  end
end

X = Q*X*Q2';


% Cholesky factor
%[Vx,Dx] = eig((X+X')/2);
%X = Vx*diag(sqrt(diag(Dx).*indx));
%indx=(diag(Dx)>0);
%X = Ech'\(Vx(:,indx)*sqrt(Dx(indx,indx)));



