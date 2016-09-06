function [runtime] = test_multiply(limit)
    n = 1;
    runtime = zeros(6,limit);
    for i=1:limit
        n = n * 2;
        A = rand(n);
        runtime(1,i) = test_time(@matrix_multiply1,A);
        runtime(2,i) = test_time(@matrix_multiply2,A);
        runtime(3,i) = test_time(@matrix_multiply3,A);
        runtime(4,i) = test_time(@matrix_multiply4,A);
        runtime(5,i) = test_time(@matrix_multiply5,A);
        runtime(6,i) = test_time(@matrix_multiply6,A);
    end
end

function [t] = test_time(func, A)
    st = cputime;
    func(A,A);
    t = cputime - st;
end

function [ C ] = matrix_multiply1(A, B)
    [~,n] = size(A);
    C = zeros(n);
    for i=1:n
        for j=1:n
            for k=1:n
				C(i,j) = C(i,j) + A(i,k)*B(k,j);
            end
        end
    end
end

function [ C ] = matrix_multiply2(A, B)
    [~,n] = size(A);
    C = zeros(n);
    for i=1:n
        for k=1:n
            for j=1:n
				C(i,j) = C(i,j) + A(i,k)*B(k,j);
            end
        end
    end
end

function [ C ] = matrix_multiply3(A, B)
    [~,n] = size(A);
    C = zeros(n);
    for j=1:n
        for i=1:n
            for k=1:n
				C(i,j) = C(i,j) + A(i,k)*B(k,j);
            end
        end
    end
end

function [ C ] = matrix_multiply4(A, B)
    [~,n] = size(A);
    C = zeros(n);
    for j=1:n
        for k=1:n
            for i=1:n
				C(i,j) = C(i,j) + A(i,k)*B(k,j);
            end
        end
    end
end

function [ C ] = matrix_multiply5(A, B)
    [~,n] = size(A);
    C = zeros(n);
    for k=1:n
        for i=1:n
            for j=1:n
				C(i,j) = C(i,j) + A(i,k)*B(k,j);
            end
        end
    end
end

function [ C ] = matrix_multiply6(A, B)
    [~,n] = size(A);
    C = zeros(n);
    for k=1:n
        for j=1:n
            for i=1:n
				C(i,j) = C(i,j) + A(i,k)*B(k,j);
            end
        end
    end
end
