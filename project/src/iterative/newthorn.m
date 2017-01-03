function [root] = newthorn(A,tol,x0,nmax,iref)
    apoly=A;
    n = length(A);
    for i=1:n-1, it=1; xn(it,i) = x0 + sqrt(-1)*x0; err=tol+1; Ndeg=n-i+1;
        if Ndeg == 1
            it=it+1; 
            xn(it,i) = -A(2)/A(1);
        else
            while it < nmax & err > tol
                [px, B] = horner(A,Ndeg,xn(it,i)); 
                [pdx, C] = horner(B,Ndeg-1,xn(it,i));
                it=it+1;
                if pdx ~= 0
                    xn(it,i) = xn(it-1,i) - px/pdx;
                    err = max(abs(xn(it,i)-xn(it-1,i)),abs(px));
                else
                    fprintf('Stop due to a vanishing p');
                    err = 0; 
                    xn(it,i) = xn(it-1,i);
                end
            end
        end
        A = B;
        if iref==1
            alfa=xn(it,i); itr=1; err=tol+1;
            while err>tol*1e-3 & itr<nmax
                [px, B] = horner(apoly,n,alfa); 
                [pdx, C] = horner(B,n-1,alfa);
                itr=itr + 1;
                if pdx ~= 0
                    alfa2 = alfa-px/pdx;
                    err = max(abs(alfa2-alfa),abs(px));
                    alfa = alfa2;
                else
                    fprintf('Stop due to a vanishing p');
                    err=0;
                end
            end
            itrefin(i)=itr-1; 
            xn(it,i)=alfa;
        end
        iter(i)=it-1;
        root(i)=xn(it,i); 
        x0 = root(i);
    end
end
