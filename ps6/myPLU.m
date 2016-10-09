function [ p, L, U ] = myPLU( A, flag )
[n,~] = size(A);
[A, p] = myLU(A, flag);
L = eye(n);
U = zeros(n);
for i=1:n
    L(i,1:i-1) = A(p(i),1:i-1);
    U(i,i:n) = A(p(i),i:n);
end
end
