function [ mu ] = wilkinson( A, n )
    a = A(n-1,n-1);
    b = A(n-1,n);
    c = A(n,n-1);
    d = A(n,n);
    u = (a+d)/2;
    v = sqrt(a*a + d*d - 2*a*d + 4*b*c)/2;
    lambda1 = u + v;
    lambda2 = u - v;
    if norm(lambda1 - A(n,n)) < norm(lambda2 - A(n,n))
        mu = lambda1;
    else
        mu = lambda2;
    end
end