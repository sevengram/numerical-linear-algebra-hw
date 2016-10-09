data = csvread('auto.csv');

b = data(:,1);
a = data(:,2);
[m,~] = size(a);
x = 0.1:0.1:250;

hold on;

% Linear
A = ones(m,2);
A(:,2) = a;
t = myLSsolve(A,b);
y = t(1) + t(2) * x;
plot(x,y,'LineWidth',1.5);
E1 = sumsqr(b-A*t)

% Quadratic
A = ones(m,3);
A(:,2) = a;
A(:,3) = a.^2;
t = myLSsolve(A,b);
y = t(1) + t(2) * x + t(3) * x.^2;
plot(x,y,'LineWidth',1.5);
E2 = sumsqr(b-A*t)

% Degree-5
A = ones(m,6);
for i = 1:5
    A(:,i+1) = a.^i;
end
t = myLSsolve(A,b);
y = t(1)+t(2)*x+t(3)*x.^2+t(4)*x.^3+t(5)*x.^4+t(6)*x.^5;
plot(x,y,'g','LineWidth',1.5);
E3 = sumsqr(b-A*t)

scatter(a,b);
axis([0, 250, 0, 60]);
legend('linear','quadratic','degree-5');
xlabel('horsepower');
ylabel('mpg');

hold off;
