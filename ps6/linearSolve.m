function [ x ] = linearSolve( A, b )
[n,~] = size(A);
[A, p] = myLU(A,1);
for i=1:n
    for j=1:i-1
        b(p(i)) = b(p(i)) - A(p(i),j) * b(p(j));
    end
end
for i=n:-1:1
    for j=i+1:n
        b(p(i)) = b(p(i)) - A(p(i),j) * b(p(j));
    end
    b(p(i)) = b(p(i)) / A(p(i), i);
end
x = b(p,:);
end
