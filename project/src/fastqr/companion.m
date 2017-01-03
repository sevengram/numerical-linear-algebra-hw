function [ U, x ] = companion( alphas )
    [n, ~] = size(alphas);
    U = gallery('tridiag',n,1,0,0);
    U(1,n) = 1;
    x = alphas;
    x(1) = x(1) + 1;
end
