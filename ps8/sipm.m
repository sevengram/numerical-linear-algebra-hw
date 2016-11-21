function [ew, ev2] = sipm(A, mu)
[n,~] = size(A);
ev1 = randn(n,1);
ev1 = ev1 ./ norm(ev1);
while 1
    ev2 = (A - mu * eye(n))\ev1;
    ev2 = ev2 ./ norm(ev2);
    ew = ev2'*A*ev2;
    if ev2(1) < 0
        ev2 = -ev2;
    end
    if norm(ev2 - ev1) < 1e-9
        break
    end
    ev1 = ev2;
end
end