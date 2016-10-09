function [ A ] = myQR( A )
    [m,n] = size(A);
    for k = 1:n;
        x = A(k:m,k);
        if x(1) == 0
            s = 1;
        else
            s = sign(x(1));
        end
        u = s*norm(x)*eye(size(x)) + x;
        u = u / u(1);
        gam = 2 / sumsqr(u);
        A(k:m,k:n) = A(k:m,k:n) - (gam*u)*(u'*A(k:m,k:n));
        A(k+1:m,k) = u(2:m-k+1);
    end
end