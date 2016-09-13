function [ b ] = myQb( A, b )
    [m,n] = size(A);
    for k=1:n
        u = eye(m-k+1,1);
        u(2:m-k+1) = A(k+1:m,k);
        gam = 2/norm(u)^2;
        b(k:m) = b(k:m) - gam*u*(u'*b(k:m));
    end
    b = b(1:n);
end
