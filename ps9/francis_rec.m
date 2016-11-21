function [ A ] = francis_rec( A, l, r )
    % l - left boundary of the sub-problem
    % r - right boundary of the sub-problem
    if l >= r
        % Only one element or empty
        return
    end
    count = 0;
    while 1
        count = count + 1;
        % Use wilkinson's shift
        mu = wilkinson(A, r);
        for i = l:r-1
            if i == l
                % Perform the shift
                x1 = A(l,l) - mu;
                x2 = A(l+1,l);
            else
                % A(i+1,i-1) is the bulge
                x1 = A(i,i-1);
                x2 = A(i+1,i-1);
            end
            c = x1/sqrt(x1^2 + x2^2);
            s = x2/sqrt(x1^2 + x2^2);
            % Left multiplication, mess the i-th and i+1-th columns
            for j = l:r
                a =  c*A(i,j) + s*A(i+1,j);
                b = -s*A(i,j) + c*A(i+1,j);
                A(i,j) = a;
                A(i+1,j) = b;
            end
            % Right multiplication, mess the i-th and i+1-th rows
            for j = l:r
                a =  c*A(j,i) + s*A(j,i+1);
                b = -s*A(j,i) + c*A(j,i+1);
                A(j,i) = a;
                A(j,i+1) = b;
            end
        end
        for i = r-1:-1:l
            % Check subdiagonal elements
            if abs(A(i+1,i)) <= 1e-16 * (abs(A(i,i)) + abs(A(i+1,i+1)))
                % Deflate
                B = francis_rec(A, l, i);
                A(l:i,l:i) = B(l:i,l:i);
                B = francis_rec(A, i+2, r);
                A(i+2:r,i+2:r) = B(i+2:r,i+2:r);
                return
            end
        end
    end
end

