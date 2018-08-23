    function p=RecursiveSearch(p,P,Psteps,j,jmax,activeDims)
        dim=activeDims(j);
        while p(dim)<P(2,dim)
            p(dim)=p(dim)+Psteps(dim);
            if j<jmax
                p=RecursiveSearch(p,P,Psteps,j+1,jmax,activeDims);
            else
                disp(p);
            end;
        end;
            p(dim)=P(1,dim);
 
%             
% function p=testRecursion(p,p0,pMax,ind,indMax)
%         while p(ind)<pMax(ind)
%             p(ind)=p(ind)+1;
%             if ind<indMax
%                 p=testRecursion(p,p0,pMax,ind+1,indMax);
%             else
%                 disp(p);
%             end;
%         end;
%         p(ind)=p0(ind);        