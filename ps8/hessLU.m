function [ A, p ] = hessLU( A, pivoting )
[n,~] = size(A);
p = 1:n;
for k=1:n-1
    if pivoting
        if A(p(k+1),k) > A(p(k), k)
            t = p(k+1);
            p(k+1) = p(k);
            p(k) = t;
        end
    end
    A(p(k+1),k) = A(p(k+1),k)/A(p(k),k);
    for j=k+1:n
        A(p(k+1),j) = A(p(k+1),j) - A(p(k+1),k) * A(p(k),j);
    end
end
end