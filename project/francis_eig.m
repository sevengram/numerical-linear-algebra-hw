function [ result ] = francis_eig( A, l, r )
    % l - left boundary of the sub-problem
    % r - right boundary of the sub-problem
    if l == r
        result = A(l,l);
        return
    end
    if l == r - 1
        [rho1, rho2] = wilkinson(A, r);
        result = [rho1, rho2];
        return
    end
    while 1
        % Use wilkinson's shift
        [rho1, rho2] = wilkinson(A, r);
        for k = l:r-1
            if k == l
                % Perform the shift
                x = zeros(r+1-l, 1);
                x(1) = rho1;
                x = (A(l:r, l:r) - rho2*eye(r+1-l))*(A(l:r, l) - x);
            else
                x = A(k:r, k-1);
            end
            if x(1) == 0
                s = 1;
            else
                s = sign(x(1));
            end
            u = s*norm(x)*eye(size(x)) + x;
            gam = 2 / sumsqr(u);
            if k == 1
                A(k:r, k:r) = A(k:r, k:r) - gam*u*(u'*A(k:r, k:r));
            else
                A(k:r, k-1:r) = A(k:r, k-1:r) - gam*u*(u'*A(k:r, k-1:r));
            end
            A(l:r, k:r) = A(l:r, k:r) - gam*(A(1:r, k:r)*u)*u';
        end
        for k = r-1:-1:l
            % Check subdiagonal elements
            if abs(A(k+1,k)) <= 1e-16 * (abs(A(k,k)) + abs(A(k+1,k+1)))
                % Deflate
                result = zeros(r-l+1,1);
                result(1:k-l+1) = francis_eig(A, l, k);
                result(k-l+2:r-l+1) = francis_eig(A, k+1, r);
                return
            end
        end
    end
end

