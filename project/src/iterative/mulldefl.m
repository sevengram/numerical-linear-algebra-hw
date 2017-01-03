function root = mulldefl(A,tol,x0,x1,x2,nmax,iref)
apoly=A;
n = length(A);
for i=1:n-1
    xn(1,i)=x0; xn(2,i)=x1; xn(3,i)=x2;
    it=1; err=tol+1; k=2; Ndeg=n-i+1;
    if Ndeg==1
        it=it+1; k=0; 
        xn(it,i)=-A(2)/A(1);
    else
        while err>tol
            k=k+1; it=it+1;
            [f0,B]=horner(A,Ndeg,xn(k-2,i)); 
            [f1,B]=horner(A,Ndeg,xn(k-1,i));
            [f2,B]=horner(A,Ndeg,xn(k,i));
            f01=(f1-f0)/(xn(k-1,i)-xn(k-2,i)); f12=(f2-f1)/(xn(k,i)-xn(k-1,i));
            f012=(f12-f01)/(xn(k,i)-xn(k-2,i));
            w=f12+(xn(k,i)-xn(k-1,i))*f012;
            arg=power(w,2)-4*f2*f012; d1=w-sqrt(arg);
            d2=w+sqrt(arg); den=max(d1,d2);
            if den ~= 0
                xn(k+1,i)=xn(k,i)-(2*f2)/den;
                err=abs(xn(k+1,i)-xn(k,i));
            else
                fprintf('Vanishing denominator');
                return
            end
        end
    end
    A=B;
    radix=xn(k+1,i);
    if iref==1
        alfa=radix; itr=1; err=tol+1;
        while err>tol
            [px,B]=horner(apoly,n,alfa); [pdx,C]=horner(B,n-1,alfa);
            if pdx == 0
                fprintf(' Vanishing derivative '); err=0;
            end
            itr=itr+1;
            if pdx ~= 0
                alfa2=alfa-px/pdx; err=abs(alfa2-alfa); alfa=alfa2;
            end
        end
        itrefin(i)=itr-1; xn(k+1,i)=alfa; radix=alfa;
    end
    iter(i)=it; root(i)=radix; 
    if Ndeg > 1
        [px,B]=horner(A,Ndeg-1,xn(k+1,i)); 
    end
end