c1 = zeros(1,11);
c2 = zeros(1,11);
for i=2:12
    A = hilb(i);
    c1(1,i-1) = norm(A,1)*norm(inv(A),1);
    c2(1,i-1) = estimate_cond_no(A);
end
x = (2:12);
figure;
semilogy(x, c1,'r', x, c2, 'b', 'LineWidth', 2);
legend('exact', 'estimate');
xlabel('n');
ylabel('condition number');