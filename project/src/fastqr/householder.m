function [ P ] = householder( x )
    u = norm(x)*[1;0;0] - x;
    gam = 2 / sumsqr(u);
    P = eye(3) - gam*(u*u');
end
