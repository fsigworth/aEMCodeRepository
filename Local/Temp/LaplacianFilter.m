function out=LaplacianFilter(in,mask,fc,iters)
n=size(in);
r=zeros(n,'single');
for i=1:iters
    out=GaussFiltDCT(r,fc);
    r=out+mask.*(in-r);
%     if mod(i,100)==0
% %          plot(sect(r));
%         imacs(r);
%         title(i);
%         drawnow;
%     end;
end;

