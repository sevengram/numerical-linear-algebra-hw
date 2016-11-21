function [ lambda1, lambda2 ] = wilkinson( A, n )
    a = A(n-1,n-1);
    b = A(n-1,n);
    c = A(n,n-1);
    d = A(n,n);
    u = (a+d)/2;
    v = sqrt(a*a + d*d - 2*a*d + 4*b*c)/2;
    lambda1 = u + v;
    lambda2 = u - v;
end