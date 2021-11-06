function [ X,res ] = lyap_solve( A,B )

%function X = lyap2solve(A,B);
% 
% Solve  A X + X A' + B = 0
%
% Bartels Stewart Method

n = size(A,1);

[Q,R]=schur(A,'complex');

C = Q'*B*Q;

X = zeros(n,n);

for i = n:-1:1
    X(i,i) = -C(i,i)/(R(i,i)+R(i,i)');
    X(1:i-1,i) = (R(1:i-1,1:i-1)+R(i,i)'*eye(i-1,i-1))\(-C(1:i-1,i)-X(i,i)*R(1:i-1,i));
    X(i,1:i-1) =  X(1:i-1,i)';
    C(1:i-1,1:i-1) = C(1:i-1,1:i-1)+X(1:i-1,i)*R(1:i-1,i)'+R(1:i-1,i)*X(1:i-1,i)';
end

X = Q*X*Q';

res = norm(A*X+X*A'+B,'fro');

end

