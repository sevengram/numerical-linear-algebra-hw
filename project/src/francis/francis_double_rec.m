function [ A ] = francis_double_rec( A, l, r )
    % l - left boundary of the sub-problem
    % r - right boundary of the sub-problem
    if l >= r - 1
        % 1x1 or 2x2 or empty
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
            A(l:r, k:r) = A(l:r, k:r) - gam*(A(l:r, k:r)*u)*u';
        end
        for k = r-1:-1:l
            % Check subdiagonal elements
            if abs(A(k+1,k)) <= 1e-16 * (abs(A(k,k)) + abs(A(k+1,k+1)))
                % Deflate
                B = francis_double_rec(A, l, k);
                A(l:k,l:k) = B(l:k,l:k);
                B = francis_double_rec(A, k+2, r);
                A(k+2:r,k+2:r) = B(k+2:r,k+2:r);
                return
            end
        end
    end
end

