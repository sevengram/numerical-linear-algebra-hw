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
        gamma = 2 / norm(u)^2;
        A(k:m,k:n) = A(k:m,k:n) - (gamma*u)*(u'*A(k:m,k:n));
        A(k+1:m,k) = u(2:m-k+1);
    end
end