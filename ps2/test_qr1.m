limit = 8;
gs_delta = zeros(1,limit);
mgs_delta = zeros(1,limit);
qr_delta = zeros(1,limit);
n = 1;
for i=1:limit
    n = n*2;
    A = rand(2*n, n);
    I = eye(n);
    [Q,~] = gram_schmidt(A);
    gs_delta(1,i) = norm(Q'*Q - I, inf);
    [Q,~] = modified_gram_schmidt(A);
    mgs_delta(1,i) = norm(Q'*Q - I, inf);
    [Q,~] = qr(A,0);
    qr_delta(1,i) = norm(Q'*Q - I, inf);
end
x_plot = 2.^(1:limit);
figure;
loglog(x_plot,gs_delta, 'b-', x_plot, mgs_delta, 'r-', x_plot, qr_delta, 'g-', 'LineWidth', 2);
legend('gram schmidt', 'modified gram schmidt', 'matlab qr','Location','southeast');
xlabel('n');
ylabel('delta');
