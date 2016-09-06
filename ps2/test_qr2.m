gs_delta = zeros(1,9);
mgs_delta = zeros(1,9);
qr_delta = zeros(1,9);
for i=2:10
    eps = 10.^(-i);
    A = [1,1,1;eps,0,0;0,eps,0;0,0,eps];
    I = eye(3);
    [Q,~] = gram_schmidt(A);
    gs_delta(1,i-1) = norm(Q'*Q - I, inf);
    [Q,~] = modified_gram_schmidt(A);
    mgs_delta(1,i-1) = norm(Q'*Q - I, inf);
    [Q,~] = qr(A,0);
    qr_delta(1,i-1) = norm(Q'*Q - I, inf);
end
x_plot = 10.^(-2:-1:-10);
figure;
loglog(x_plot,gs_delta, 'b-', x_plot, mgs_delta, 'r-', x_plot, qr_delta, 'g-', 'LineWidth', 2);
legend('gram schmidt', 'modified gram schmidt', 'matlab qr','Location','southeast');
set(gca, 'xdir', 'reverse');
xlabel('epsilon');
ylabel('delta');
