function test_matrix(limit)
    n = 1;
    time_col = zeros(1,limit);
    time_row = zeros(1,limit);
    for i=1:limit
        n = n * 2;
        x = rand(n, 1);
        A = rand(n);
        tic;
        axrow(A,x);
        time_row(1,i) = toc;
        tic;
        axcol(A,x);
        time_col(1,i) = toc;
    end
    x_plot = 2.^(1:limit);
    figure;
    plot(x_plot,time_col, 'b--', x_plot, time_row, 'r-','LineWidth', 2);
    legend('axcol','axrow');
    xlabel('n');
    ylabel('Runtime (seconds)');
end

function y = axcol(A, x)
    [n,~] = size(x);
    y = zeros(n,1);
    for j=1:n
        for i=1:n
            y(i,1) = y(i,1) + A(i,j)*x(j,1);
        end
    end
end

function y = axrow(A, x)
    [n,~] = size(x);
    y = zeros(n,1);
    for i=1:n
        for j=1:n
            y(i,1) = y(i,1) + A(i,j)*x(j,1);
        end
    end
end