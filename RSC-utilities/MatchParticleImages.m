% MatchParticleImages.m
% Given two si files, read the particle images and match them by
% correlation.  Then update the active flags of the second si file, anding
% them with the corresponding active flags from the first.
[siFiles, pa]=uigetfile({'*si.mat' '*si.txt'},'Select two si files','multiselect','on');
if isnumeric(pa)
    return
end;
cd(pa);
if ~iscell(siFiles)
    siFiles={siFiles};
end;
if numel(siFiles)<2
    disp('Need to select a pair of si files.');
    return
end;
disp('Reading:');
disp(siFiles{1});
si1=LoadStruct(siFiles{1});
disp([num2str(numel(si1.miIndex)) ' particles.']);
disp(siFiles{2});
si2=LoadStruct(siFiles{2});
disp([num2str(numel(si2.miIndex)) ' particles.']);
stFiles=cell(1,2);
for i=1:2
    [pa,nm,ex]=fileparts(siFiles{i});
    stFiles{i}=[nm(1:end-1) 'tack.mrc'];
end;
disp('Reading:');
m1=ReadMRC(stFiles{1});
disp(stFiles{1});
m2=ReadMRC(stFiles{2});
disp(stFiles{2});

[n,n,nim1]=size(m1);
nim2=size(m2,3);

m1=reshape(m1,n^2,nim1);
m2=reshape(m2,n^2,nim2);

disp('Correlating');

r=m1'*m2;  % m1 is rows, m2 is columns.

%%
figure(10);
imags(r);
xlabel([stFiles{1} '  ' num2str(nim1)],'interpreter', 'none');
ylabel([stFiles{2} '  ' num2str(nim2)],'interpreter', 'none');

%%
% For each m2 image, find the index of the corresponding m1 image if any.
% The index is inds1; is none match, the entry in inds1 is zero.
thresh=.2;
inds1=zeros(nim2,1);
vals1=zeros(nim2,1);
for i=1:nim2
    [vals1(i),inds1(i)]=max(r(:,i));
end;
inds1(vals1<thresh)=0;
% [v1Sort,j]=sort(vals1,1,'descend');
% valid=false(nim2,1);
% valid(j(1:nim1))=true;
% inds1(~valid)=0;

% For each m1 image, find the index of the corresponding m2 image if any.
inds2=zeros(nim1,1);
vals2=zeros(nim1,1);
for i=1:nim1
    [vals2(i),inds2(i)]=max(r(i,:));
    
end;
inds2(vals2<thresh)=0;

%%

disp(' ');
disp([stFiles{1} ' active flags:']);
naf1=size(si1.activeFlags,2);
for i=1:naf1
    if numel(si1.activeFlagLog)<i
        si1.activeFlagLog{i,1}='';
    end;
    disp([num2str([i sum(si1.activeFlags(:,i))]) '  ' si1.activeFlagLog{i}]);
end;
k1=naf1;
k1=MyInput('Select data at selection ',k1);
if k1<1 % reset the active flags
    si1.activeFlags=true(size(si1.miIndex));
    k1=1;
    disp('Re-initializing the active flags for si1.');
end;
active1=si1.activeFlags(:,k1);
active1c=active1;
active2=false(nim2,1);
for i=1:nim2
    if inds1(i)>0 
        active2(i)=active1c(inds1(i));
        active1c(inds1(i))=false;  % avoid duplicates
    end;
end;
disp('adding to stack2 active flags');
k2=size(si2.activeFlags,2);
si2.activeFlags(:,k2+1)=si2.activeFlags(:,k2) & active2;
disp([num2str(sum(si2.activeFlags(:,k2+1))) ' active flags for si2']);
logStr=[date '  MatchParticleImages of set ' num2str(k2)...
    ' with ' siFiles{1} ' set ' num2str(k1)];
si2.activeFlagLog{k2+1,1}=logStr;
disp(logStr);

str=input('Write the si file [n]? ','s');
if numel(str)>0 && lower(str(1))=='y'
    si=si2;
    save(siFiles{2},'si');
    disp([siFiles{2} ' written.']);
end;





return


inds1(~valid)=0;  % inds1(i) is the index of m1 which matches m2(i)
%%
% si1=load('sq10_350p128vtsi.mat');
% si2=load('sq10_350m2p128v2tsi.mat');
% %%
naf=size(si1.si.activeFlags,2);
si2.si.activeFlags=false(nim2,naf);
for i=1:nim2
    if inds1(i)>0
        si2.si.activeFlags(i,:)=si1.si.activeFlags(inds1(i),:);
    end;
end;
si2.si.activeFlagLog=si2.si.activeFlagLog(:);
si2.si.activeFlagLog{naf,1}='MatchParticleImages from_sq10_350p128vtsi.mat';

%%
si=si2.si;
save('sq10_350m2p128v2tsi.mat','si');