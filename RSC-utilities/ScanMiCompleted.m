% ScanMiCompleted.m
% Check for processed files matching mi files.  Any missing are put into a list.
%%
% infoPattern='mi.txt';
% procPattern='rscc.mat';
% pa2='Merged/';
% % 
% % 
% % Have the user select some mi files
% disp('Select the set of mi files');
% 
% [fname1, pa1]=uigetfile({'*mi.txt' '*mi.mat'},'Select mi files to examine','multiselect','on');
% if ~iscell(fname1)
%     fname1={fname1};
% end;
% [rootDir1, infoDir1]=ParsePath(pa1);
% cd(rootDir1);
% pa1=AddSlash(pa1);
% 
% nmi=numel(fname1);
% baseName1=cell(nmi,1);
% for i=1:nmi
%     p=strfind(fname1{i},infoPattern);
%     baseName1{i}=fname1{i}(1:p-1);
% end;
    
%% scan for the processed files
d=dir(pa2);
k=0;
baseName2=cell(0,1);
for i=1:numel(d)
    p=strfind(d(i).name,procPattern);
    if numel(p)>0
        k=k+1;
        baseName2{k,1}=d(i).name(1:p-1);
    end;
end;

%%
% pick up all the matches
matches=false(nmi,1);
for i=1:nmi
    p=strcmp(baseName1{i},baseName2);
    if any(p)
        matches(i)=true;
    end;
end

disp([num2str(sum(matches)) ' matches.']);

toDoInfos=fname1(~matches);
whos toDoInfos
save toDoInfos.mat toDoInfos

return

nim2=numel(fname2);
ptrs21=zeros(nim2,1);

for i=1:nim2
    p=strcmp(fname2{i},fname1);
    if any(p)
        ptrs21(i)=find(p,1);
    end;
end;    

sum(ptrs12>0)

sum(ptrs21>0)

%%  Create the common list of names
fileNames=fname1(ptrs12>0);
nim=numel(fileNames);

mis1=cell(nim,1);
for i=1:nim
    disp(fileNames{i});
    mi=ReadMiFile([pa1 fileNames{i}]);
    mis1{i,1}=mi;
end;

save([pa1 'mis1.mat'],'mis1');

mis2=cell(nim,1);
for i=1:nim
    disp(fileNames{i});
    mi=ReadMiFile([pa2 fileNames{i}]);
    mis2{i,1}=mi;
end;

save([pa2 'mis2.mat'],'mis2');

%%
% save AllMis.mat mis
return

%% Pick up defocus values
nim=numel(fileNames);
defs=zeros(nim,2);
astig=zeros(nim,2);
for i=1:nim
    defs(i,1)=mis1{i}.ctf(1).defocus;
    defs(i,2)=mis2{i}.ctf(1).defocus;
    astig(i,1)=mis1{i}.ctf(1).deltadef;
    astig(i,2)=mis2{i}.ctf(1).deltadef;
end;

%%
% Get amplitudes
ds=4;
nv=0;
for i=1:nim
    mi=mis{i,1};
    niv=numel(mi.vesicle.x);
    % multiply the vesicle model by the voxel size.
    
    if niv>0 && isfield(mi.vesicle,'extraPeaks') && numel(mi.vesicle.extraPeaks)>0
        maxR=max(mi.vesicle.r(all(mi.vesicle.ok,2),1))/ds;
        n=NextNiceNumber(maxR*3);
        vx=n/2+1;
        vy=vx;
        vd=meDownsampleVesicleModel(mi.vesicleModel,ds)*ds*mi.pixA;
        
        for ind=1:niv
            if all(mi.vesicle.ok(ind,:))
                nv=nv+1;  % count of all vesicles
                % Get the coordinates and radius, scaled down by ds
                vr=mi.vesicle.r(ind,:)/ds;
                
                % Accumulate the vesicle density
                %     sumv=sumv-mi.vesicle.s(ind)*VesicleFromModel(n,vr,vd,[vx vy]);
                v0=VesicleFromModelGeneral(n,vr,vd,[vx vy],mi.vesicle.s(ind,:,1));
                
                % -------------Create extra ring components--------------
                % extra amps are part of the s array
                exAmps=ds*mi.pixA*mi.vesicle.s(ind,:,2:end);
                exPos=mi.vesicle.extraPeaks/ds;
                exSigma=mi.vesicle.extraSD/ds;
                v1=v0+VesicleFromRings(n,exPos,exSigma,vr,[vx vy],exAmps);
                mulr=50;
                addr=0;
                subplot(221);
                imaga(v0*mulr+addr);
                subplot(222);
                imaga(v1*mulr+addr);
                subplot(224);
                imaga((v1-v0)*mulr+128);
                subplot(223);
                plot([sect(v0) sect(v1)]);
                axis([0 n 0 5]);
                drawnow;
            end;
        end;
    end;
end;



% <</Users/fred/EMWork/Hideki/151117/KvLipo80slot3/Trends/VesicleRings.jpg>>



% % Get amplitudes
% nv=0;
% for i=1:nim
%     ves=v{i,1};
%     niv=numel(ves.x);
%     for j=1:niv
%         if all(ves.ok(j,:))
%             nv=nv+1;
%     ves0=ves;
%     ves0.extraPeaks=[];
%     v=meMakeModelVesicles(mi,n,j,0,0);
%     v=
%     s1(nv)=ves.s(
% e1=zeros(nim,nbins);
% eAll=zeros(nim,nbins);
% figure(1);
% subplot(111);
% for i=1:nim
%     %     load(fname{i});
%     mi=ReadMiFile(fname{i});
%     if numel(mi.vesicle.x)>0
%             okVes=all(mi.vesicle.ok,2);
%             radii=mi.vesicle.r(:,1)*mi.pixA;
% %             radii=mi.vesicle.r(:,1);
%             val(i,:)=hist(radii,xbins);
%
%     bar(xbins,val(i,:));
%     title(fname{i},'interpreter','none');
%     drawnow;
%     end;
% end;
%
% bar(xbins,sum(val,1));
% xlabel('Vesicle radius, Å');
% ylabel('Frequency');
