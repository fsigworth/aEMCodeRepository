% xShowClassEvolution.m
% Given a reconstruction directory, read the moi files and
% plot the moi.pVols vs iteration number.
% 

figure(20);
d=dir;
i=0;
as=zeros(0,3);
for ind=3:numel(d)
        a=sscanf(d(ind).name,'%c%d%s'); % e.g. i01a.mrc
        if numel(a)==11 && a(5)=='m' % got all fields
            i=i+1;
            as(i,:)=[a(2) a(3)-96 ind];
        end;
end;
firstIter=min(as(:,1))
lastIter=max(as(:,1))
load(d(as(1,3)).name);
pVols=moi.pVols';
nVols=numel(pVols);
iters=as(1,1);
if iters==1
    pVols(1,:)=NaN;   % don't plot the first point.
end;
for i=2:size(as,1)
    name=d(as(i,3)).name;
    disp(name);
    load(name);
    pVols(i,:)=moi.pVols';
    iters(i,1)=as(i,1);
end;
%%
plot(iters,pVols);
legend(num2str((1:nVols)'));
xlabel('Iteration');
ylabel('Fraction of particles');

return
%%









% firstIter=lastIter;

for ind=firstIter:lastIter
    is=find(as(:,1)==ind);
    if numel(i)<1 break; end;
    %%
    for i=is'
        figure(i-is(1)+1);
        vi=i-is(1)+1;
        name=d(as(i,4)).name;
        disp(name);
        if vi==1
            v=ReadMRC(name);
        else
        v(:,:,:,vi)=ReadMRC(name);
        end;
        ShowSections(v(:,:,:,vi),[],45);
        title(name);
    end;
    pause;
end;


return

       
iTxt=['i' sprintf('%02d',ind)];
 [v1,s]=ReadMRC(['mrc/' iTxt 'av01.mrc']);
 v2=ReadMRC(['mrc/' iTxt 'bv01.mrc']);
 [vol2ali,tz,gamma,mirror]=reAlignVolumes(v1,v2);
 vol2ali=v2;
 
disp(['iter ' sprintf('%3d',ind) '  alignment: ' sprintf('[%5.1f %5.1f  %1d ]',[tz gamma mirror])]);

 n=size(v2,1);
 cn=ceil(n/2+1);
%   vmsk=fuzzymask(n,3,0.2*n,0.1*n,[cn cn cn-12]);
%   vmsk=fuzzymask(n,3,0.3*n,0.1*n,[cn cn cn]);
vmsk=min(1,fuzzymask(n,3,[.32 .32 .15]*n,.1*n,[cn cn cn+.18*n])+fuzzymask(n,3,0.22*n,0.1*n,[cn cn cn-.06*n])); figure(1);
 ShowSections(vmsk.*(v1),[],45);
 [p2,p1]=ParsePath(pwd);
 n=size(v1,1);
 fs=(0:n/2-1)/(n*s.pixA);
 fsc=FSCorr(vmsk.*v1,vmsk.*vol2ali);
 subplot(339);
 plot(fs,[fsc fsc*0]);
 xlabel('Frequency, Å^{-1}');
 ylabel('FSC');
 title([p1 '  ' iTxt]);
 
 figure(2);
 ShowSections(vmsk.*vol2ali,[],45);
 
 pause
return


load ri
 subplot(332);
 title(ri.siName,'interpreter','none');
 
index=3;
fscs(:,index)=fsc;
