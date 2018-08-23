gammas=0:11.25:90;
angles(:,3)=gammas';
angles(:,2)=90;
projs=rsMakeTemplates(angles,map);
%%
nang=size(angles,1);
while 1
for i=1:nang
    imacs(projs(:,:,i));
    title(angles(i,3));
    pause
end;
end;