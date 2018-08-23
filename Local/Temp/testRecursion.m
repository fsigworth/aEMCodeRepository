function p=testRecursion(p,p0,pMax,ind,indMax)
        while p(ind)<pMax(ind)
            p(ind)=p(ind)+1;
            if ind<indMax
                p=testRecursion(p,p0,pMax,ind+1,indMax);
            else
                disp(p);
            end;
        end;
        p(ind)=p0(ind);        