function [ result ] = francis_roots( alphas )
    [~, n] = size(alphas);
    alphas = alphas / alphas(1);
    A = gallery('tridiag',n-1,1,0,0);
    A(:,n-1) = -alphas(n:-1:2)';
    result = francis_eig(A,1,n-1);
end
