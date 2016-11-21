function [ B ] = wrap( A, k, n )
    [s,~] = size(A);
    B = blkdiag(eye(k-1), A, eye(n-k-s+1));
end
