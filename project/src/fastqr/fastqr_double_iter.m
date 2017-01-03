function [ V, W, x, y] = fastqr_double_iter( V, W, x, y, mu1, mu2)
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
    
    % Step 2
    beta = x(n);
    for k = n-1:-1:4
        beta = V(:,1,k)'*[x(k); beta];
    end
    
    VV = eye(n);
    WW = eye(n);
    for i = n-1:-1:2
        VV = VV * wrap(V(:,:,i), i, n);
    end
    for i = 1:n-2
        WW = WW * wrap(W(:,:,i), i, n);
    end
    A = VV*WW - x*y'
    
    % Step 3
    P = zeros(3,3,n);
    rho = zeros(n,1);
    delta = zeros(n,1);
    theta = 0;
    xi = 0;
    eta = 0;
    for k = 1:n-2
        if k == 1
            A1 = [g(:,1)' * h(:,1); sigma(1); x(3:n) * y(1)] - x(1:n,1) * y(1);
            A13 = zeros(3,n);
            for j = 1:3
                HH = zeros(2,n);
                BB = eye(2);
                HH(:,j) = h(:,j);
                for i = j+1:n
                    BB = BB * B(:,:,i);
                    HH(:,i) = BB * h(:,i);
                end
                A13(j,:) = g(:,j)' * HH - y(1:n,1)' * x(j);
            end
            A13(:,1) = A1(1:3);
            A13(3,2) = sigma(2) - x(3)*y(2);
            A13(1,1) = A13(1,1) - mu1;
            A13(2,2) = A13(2,2) - mu1;
            A13(3,3) = A13(3,3) - mu1;
            A1(1) = A1(1) - mu2;
            
            T = [A13 * A1, [g(:,1)' * h(:,1) - x(1)*y(1); sigma(1) - x(2)*y(1); 0]];
            R = [g(:,1)'*h(:,1), g(:,1)'*B(:,:,2)*h(:,2), g(:,1)'*B(:,:,2)*B(:,:,3)*h(:,3), g(:,1)'*B(:,:,2)*B(:,:,3)*B(:,:,4)*h(:,4);
                 sigma(1),       g(:,2)'*h(:,2),          g(:,2)'*B(:,:,3)*h(:,3),          g(:,2)'*B(:,:,3)*B(:,:,4)*h(:,4);
                 x(3)*y(1),      sigma(2),                g(:,3)'*h(:,3),                   g(:,3)'*B(:,:,4)*h(:,4);
                 x(4)*y(1),      x(4)*y(2),               sigma(3),                         g(:,4)'*h(:,4)];
        elseif 1 < k && k < n-2
            T = [sigma(k-1), rho(k); theta, sigma(k); xi, eta] - x(k:k+2) * y(k-1:k)';
            R = [rho(k),         delta(k),                gbh1,                             g(:,k)'*B(:,:,k+2)*B(:,:,k+3)*h(:,k+3);
                 sigma(k),       rho(k+1),                gbh2,                             g(:,k+1)'*B(:,:,k+2)*B(:,:,k+3)*h(:,k+3);
                 eta,            sigma(k+1),              g(:,k+2)'*h(:,k+2),               g(:,k+2)'*B(:,:,k+3)*h(:,k+3);
                 x(k+3)*y(k),    x(k+3)*y(k+1),           sigma(k+2),                       g(:,k+3)'*h(:,k+3)];
        else
            T = [sigma(k-1), rho(k); theta, sigma(k); xi, eta] - x(k:k+2) * y(k-1:k)';
            R = [rho(k),         delta(k),                gbh1;
                 sigma(k),       rho(k+1),                gbh2;
                 eta,            sigma(k+1),              g(:,k+2)'*h(:,k+2)]; 
        end
        P(:,:,k) = householder_givens(T);
        if k < n-2
            T = blkdiag(P(:,:,k)',1) * R * blkdiag(P(:,:,k),1);
        else
            T = P(:,:,k)' * R * P(:,:,k);
        end
        sigma(k) = T(2,1);
        rho(k+1) = T(2,2);
        delta(k+1) = T(2,3);
        theta = T(3,1);
        sigma(k+1) = T(3,2);
        rho(k+2) = T(3,3);
        if k < n-2
            gbh1 = T(2,4)
            bh = B(:,:,k+3) * h(:,k+3);
            g(:,k+1) = [1;(gbh1 - bh(1))/bh(2)];
            gbh2 = T(3,4)
            g(:,k+2) = [1;(gbh2 - bh(1))/bh(2)]
            xi = T(4,1);
            eta = T(4,2);
            sigma(k+2) = T(4,3);
        end
        x(k:k+2) = P(:,:,k)'*x(k:k+2);
        y(k:k+2) = P(:,:,k)'*y(k:k+2);
        if k == 1
            x_2 = x(2);
            x_3 = x(3);
        end
    end
    
    P(2:3,2:3,n-1) = givens(sigma(n-2) - x(n-1)*y(n-2), theta - x(n)*y(n-2));
    P(1,1,n-1) = 1;
    x(n-1:n) = P(2:3,2:3,n-1)'*x(n-1:n);
    y(n-1:n) = P(2:3,2:3,n-1)'*y(n-1:n);
    
    % Step 4
    [V_3, beta] = givens(x_3, beta);
    [V_2, ~] = givens(x_2, beta);
    D = blkdiag(1,1,V_3')*blkdiag(1,V_2',1)*blkdiag(P(:,:,1)',1)*blkdiag(1,1,V(:,:,3))*blkdiag(1,V(:,:,2),1)*blkdiag(W(:,:,1),1)*blkdiag(1,W(:,:,2));
    D
    V(:,:,2) = V_2;
    V(:,:,3) = V_3;

    
end
