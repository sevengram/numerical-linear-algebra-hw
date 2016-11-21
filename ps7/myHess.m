function [ A ] = myHess( A )
    [~,n] = size(A);
    for k = 1:n-2;
        x = A(k+1:n,k);
        if x(1) == 0
            s = 1;
        else
            s = sign(x(1));
        end
        u = s*norm(x)*eye(size(x)) + x;
        gam = 2 / sumsqr(u);
        A(k+1:n,k:n) = A(k+1:n,k:n) - gam*u*(u'*A(k+1:n,k:n));
        A(1:n, k+1:n) = A(1:n, k+1:n) - gam*(A(1:n, k+1:n)*u)*u';
    end
end

