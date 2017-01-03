function [ V, W, x, y] = fastqr_single_iter( V, W, x, y, mu)
    [n, ~] = size(x);
    
    % Step 1
    % Quasiseparable parametrization:  h, g, B, sigma
    h = zeros(2,n);
    B = zeros(2,2,n);
    g = zeros(2,n);
    sigma = zeros(n-1,1);
    
    h(:,1:n-2) = W(1:2,1,1:n-2);
    B(:,:,2:n-1) = W(1:2,2:3,1:n-2);
    h(:,n-1) = [1;0];
    h(:,n) = [0;1];
    B(:,:,n) = eye(2);
    
    g(:,1) = [1;0];
    rr = [0,1];
    for k = 1:n-2
        T = V(:,:,k+1)*[rr,0;0,0,1]*W(:,:,k);
        sigma(k) = T(1,1);
        g(:,k+1) = T(1,2:3)';
        rr = T(2,2:3);
    end
    sigma(n-1) = rr * h(:,n-1);
    g(:,n)= rr';

    VV = eye(n);
    WW = eye(n);
    for i = n-1:-1:2
        VV = VV * wrap(V(:,:,i), i, n);
    end
    for i = 1:n-2
        WW = WW * wrap(W(:,:,i), i, n);
    end
    U = VV*WW
    A = U - x*y'
    eig(A)
    
    % Step 2
    beta = x(n);
    for k = n-1:-1:3
        beta = V(:,1,k)'*[x(k); beta];
    end
    
    % Step 3
    r = zeros(n,1);
    d = zeros(n,1);
    P = zeros(2,2,n);
    for k = 1:n-1
        if k == 1
            [ P(:,:,k), ~ ] = givens(g(:,1)' * h(:,1) - x(1)*y(1) - mu, sigma(1) - x(2)*y(1));
            sigma(1) =  P(1,2,k) * g(:,1)' * h(:,1)   + P(2,2,k) * sigma(1);
        else
            [ P(:,:,k), ~ ] = givens(r(k)             - x(k)*y(k-1)   , r(k+1)   - x(k+1)*y(k-1));
            sigma(k) =  P(1,2,k) * d(k)               + P(2,2,k) * d(k+1);
        end
        g(:,k+1) = (P(1,2,k) * g(:,k)' * B(:,:,k+1) + P(2,2,k) * g(:,k+1)')';
        if k ~= n-1
            T = [sigma(k) , g(:,k+1)' * h(:,k+1) ; x(k+2)*y(k) , sigma(k+1)] * P(:,:,k);
            r(k+1:k+2) = T(1:2,1);
            d(k+1:k+2) = T(1:2,2);
        end
        x(k:k+1) = P(:,:,k)'*x(k:k+1);
        y(k:k+1) = P(:,:,k)'*y(k:k+1);
        if k == 1
            x_2 = x(2);
        end
    end
    
    % Step 4
    [V_2, ~] = givens(x_2, beta);
    W(:,:,1) = blkdiag(1, V_2') * blkdiag(P(:,:,1)', 1) * blkdiag(1, V(:,:,2)) * W(:,:,1);
    V(:,:,2) = V_2;
    for k = 1:n-3
        X = blkdiag(P(:,:,k+1)', 1) * blkdiag(1, V(:,:,k+2)) * blkdiag(V(:,:,k+1), 1);
        Y = blkdiag(W(:,:,k), 1) * blkdiag(1, W(:,:,k+1)) * blkdiag(P(:,:,k), 1, 1);
        [Z, ~] = givens(X(1,2), X(1,3));
        C = X * blkdiag(1,Z);
        D = blkdiag(1,1,Z') * Y;
        
        [G,~] = givens(C(1,1),C(1,2));
        C = C * blkdiag(G,1);
        V(:,:,k+1) = G';
        V(:,:,k+2) = C(2:3,2:3);
        H = householder(D(1:3,1));
        D = blkdiag(H,1) * D;
        W(:,:,k) = H';
        W(:,:,k+1) = D(2:4,2:4);
    end
    V(:,:,n-1) = P(:,:,n-1)' * V(:,:,n-1);
    W(:,:,n-2) = W(:,:,n-2) * blkdiag(P(:,:,n-2),1) * blkdiag(1,P(:,:,n-1));
end

