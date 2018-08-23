
% rePurgeClasses.m

[si, imgs, names, siPath]=reLoadStackFiles;
cd(siPath);
if ~isfield(si,'activeFlags')
    si.activeFlags=true(numel(si.miIndex),1);
end;
%%
[clsName,clsPath]=uigetfile('*cls.mat','Select the classes file');
clsPath=AddSlash(clsPath);
cname=[clsPath clsName];
% name=names{end};
% clsName=[name(1:end-5) 'cls.mat'];
load(cname);  % loads 'partCls','clsLists','clsMeans'
n=size(clsMeans,1);
%%
% % DeNovo class mean calculation
% n=size(clsMeans,1);
% 
% nCls=max(size(clsMeans,3),max(partCls));
% clsMeans=zeros(n,n,nCls,'single');
% clsN=zeros(nCls,1);
% nim=numel(partCls);
% for i=1:nim
%     ind=partCls(i);
%     clsMeans(:,:,ind)=clsMeans(:,:,ind)+imgs(:,:,i);
%     clsN(ind)=clsN(ind)+1;
% end;
% for j=1:nCls
%     if clsN(j)>0
%         clsMeans(:,:,j)=clsMeans(:,:,j)/sqrt(clsN(j));
%     end;
% end;

dsMag=2;
ImagicDisplay3(ExpandImage(clsMeans,dsMag));
drawnow;
b=0;
selClasses=[];
ok=si.activeFlags(:,end);
disp(['Initial count: ' num2str(sum(ok))]);

while b~='q'
    [ind,loc,b]=ImagicDisplay3('GetClick');
    if numel(b)<1
        b=0;
    end;
    switch b
        case 1
        mk.marker='bs';
        mk.markerSize=n*1.4;
        ImagicDisplay3('Mark',ind,mk);
        disp(numel(clsLists{ind}));
        case 2
        mk.marker='ro';
        mk.markerSize=n*1.8;
        ImagicDisplay3('Mark',ind,mk);
        selClasses(end+1)=ind;
        ok(partCls==ind)=false;
%         ok(clsLists{ind})=false;
        disp(['Remaining: ' num2str(sum(ok))]);
    end;
    pause(0.1);
end;
disp('done');
selClasses
%%
si.activeFlags(:,end+1)=ok;
si.activeFlagLog{end+1}=[date '  rePurgeClasses'];

finalImgs=sum(ok)

str=input('Write the si file [n]? ','s');
if numel(str)>0 && lower(str(1))=='y'
    save(names{end},'si');
    disp([names{end} ' written.']);
end;

