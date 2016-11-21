function [ A ] = francis( A, mode )
    [n,~] = size(A);
    count = 0;
    while 1
        count = count+1;
        if mode == 1
            mu = 0;
        elseif mode == 2
            mu = A(n,n);
        else
            mu = wilkinson(A, n);
        end
        for i=1:n-1
            if i == 1
                % Perform the shift
                x1 = A(1,1) - mu;
                x2 = A(2,1);
            else
                x1 = A(i,i-1);
                x2 = A(i+1,i-1);
            end
            c = x1/sqrt(x1^2 + x2^2);
            s = x2/sqrt(x1^2 + x2^2);
            % Left multiplication, mess the i-th and i+1-th columns
            for j=1:n
                a =  c*A(i,j) + s*A(i+1,j);
                b = -s*A(i,j) + c*A(i+1,j);
                A(i,j) = a;
                A(i+1,j) = b;
            end
            % Right multiplication, mess the i-th and i+1-th rows
            for j=1:n
                a =  c*A(j,i) + s*A(j,i+1);
                b = -s*A(j,i) + c*A(j,i+1);
                A(j,i) = a;
                A(j,i+1) = b;
            end
        end
        if abs(A(n,n-1)) < 1e-14
            break
        end
    end
end

