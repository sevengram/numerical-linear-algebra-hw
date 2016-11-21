function [ev, ew] = sipm(A, mu, limit)
[n,~] = size(A);
ev1 = randn(n,1);
ev1 = ev1 ./ norm(ev1);
[L,U] = lu(A - mu * eye(n));
count = 0;
while 1
    ev2 = U\(L\ev1);
    ev2 = ev2 ./ norm(ev2);
    ew = ev2'*A*ev2;
    if ev2(1) < 0
        ev2 = -ev2;
    end
    count = count + 1;
    if count == limit
        break
    end
    if limit <= 0 && norm(ev2 - ev1) < 1e-16
        break
    end
    ev1 = ev2;
end
ev = ev2;
end