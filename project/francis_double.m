function [ A ] = francis_double(A)
    [n,~] = size(A);
    A = hess(A);
    [rho1, rho2] = wilkinson(A,n);
    for k=1:n-1
        if k == 1
            x = zeros(n,1);
            x(1) = rho1;
            x = (A - rho2*eye(n))*(A(:,1) - x);
        else
            x = A(k:n, k-1);
        end
        if x(1) == 0
            s = 1;
        else
            s = sign(x(1));
        end
        u = s*norm(x)*eye(size(x)) + x;
        gam = 2 / sumsqr(u);
        if k == 1
            A(k:n, k:n) = A(k:n, k:n) - gam*u*(u'*A(k:n,k:n));
        else
            A(k:n, k-1:n) = A(k:n, k-1:n) - gam*u*(u'*A(k:n,k-1:n));
        end
        A(1:n, k:n) = A(1:n, k:n) - gam*(A(1:n, k:n)*u)*u';
    end
end

