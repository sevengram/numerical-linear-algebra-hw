function [ P ] = householder_givens( B )
    H = householder(B(:,1));
    T = H' * B(:,2);
    G = givens(T(2), T(3));
    P = H * blkdiag(1,G);
end
