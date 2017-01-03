function [pnz,b] = horner(a,n,z)
b(1)=a(1);
for j=2:n
    b(j)=a(j)+b(j-1)*z;
end
pnz=b(n);
return