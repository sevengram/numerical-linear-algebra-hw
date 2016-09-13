function [ x ] = myLSsolve( A, b )
    [m,n] = size(A);
    A = myQR(A);
    c = myQb(A, b);
    R = triu(A);
    x = R(1:n,:)\c;
end
