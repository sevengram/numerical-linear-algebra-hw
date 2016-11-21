function [ output_args ] = fast_qr_iter( V, W, x, y, mu)
    % Step 1
    % Quasiseparable parametrization:  h, g, B, sigma
    h = zeros(2,n);
    B = zeros(2,2,n);
    g = zeros(2,n);
    sigma = zeros(n-1,1);
    for k = 1: n-2
        h(:,k) = W(1:2,1,k);
        B(:,:,k+1) = W(1:2,2:3,k);
    end
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
    
    % Step 3
    r = zeros(n,1);
    d = zeros(n,1);
    for k = 1:n-1
        if k == 1
            [ P, ~ ] = givens(g(:,1)' * h(:,1) - x(1)*y(1) - mu, sigma(1) - x(2)*y(1));
            sigma(1) =  P(1,2) * g(:,1)' * h(:,1)   + P(2,2) * sigma(1);
        else
            [ P, ~ ] = givens(r(k)             - x(k)*y(k-1)   , r(k+1)   - x(k+1)*y(k-1));
            sigma(k) =  P(1,2) * d(k)               + P(2,2) * d(k+1);
        end
        g(:,k+1) = (P(1,2) * g(:,k)' * B(:,:,k+1) + P(2,2) * g(:,k+1)')';
        if k ~= n-1
            T = [sigma(k) , g(:,k+1)' * h(:,k+1) ; x(k+2)*y(k) , sigma(k+1)] * P;
            r(k+1) = T(1,1);
            r(k+2) = T(2,1);
            d(k+1) = T(1,2);
            d(k+2) = T(2,2);
        end
        x(k:k+1) = P'*x(k:k+1);
        y(k:k+1) = P'*y(k:k+1);
    end
    
    % TODO: Step 4
    
end

