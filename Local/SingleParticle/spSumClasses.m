function [sums nMembers]=spSumClasses(imgs,sp)
% function [sums nMembers]=spSumClasses(imgs,sp)

[n ny nim]=size(imgs);
ncls=max(sp.class);
nMembers=zeros(ncls,1);
sums=zeros(n,ny,ncls);

for i=1:nim
    j=sp.class(i);
    nMembers(j)=nMembers(j)+1;
    sums(:,:,j)=sums(:,:,j)+imgs(:,:,i);
end;
