function [ P, b ] = givens( x1, x2 )
    c = x1/sqrt(x1^2 + x2^2);
    s = x2/sqrt(x1^2 + x2^2);
    b = c * x1 + conj(s) * x2;
    P = [c, -s; conj(s), c];
end
