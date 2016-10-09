function [ A, p ] = myLU( A, pivoting )
    [n,~] = size(A);
    p = 1:n;
    for k=1:n
        if pivoting == 1
            for i=k+1:n
                if A(p(i),k) > A(p(k), k)
                    t = p(i);
                    p(i) = p(k);
                    p(k) = t;
                end
            end
        end
        for i=k+1:n
            A(p(i),k) = A(p(i),k)/A(p(k),k);
            for j=k+1:n
                A(p(i),j) = A(p(i),j) - A(p(i),k) * A(p(k),j);
            end
        end
    end
end
