function [ V, W ] = generators( U, x )
    [n,~] = size(x);
    beta = x(n);
%   Ut = U
    
    % Decomposition of U by V and W
    V = zeros(2,2,n-1);
    W = zeros(3,3,n-2);
    for k = n-1:-1:2
        [V(:,:,k), beta] = givens(x(k), beta);
        U(k:k+1,:) = V(:,:,k)' * U(k:k+1,:);
    end
    for k=1:n-3
        t = U(k:k+2, k);
        u = - norm(t)*[1;0;0] + t;
        gam = 2 / sumsqr(u);
        W(:,:,k) = eye(3) - gam*(u*u');
        U(k:k+2, k:n) = U(k:k+2, k:n) - gam*u*(u'*U(k:k+2, k:n));
    end
    W(:,:,n-2) = U(n-2:n, n-2:n);
     
%     VV = eye(n);
%     WW = eye(n);
%     for i = n-1:-1:2
%         VV = VV * wrap(V(:,:,i), i, n);
%     end
%     Ut = VV'*Ut
%     for i = 1:n-2
%         WW = WW * wrap(W(:,:,i), i, n);
%     end
%     Ut = WW'*Ut
end

