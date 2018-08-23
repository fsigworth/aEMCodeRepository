% rlRestoreGroups
% Put the micrograph number into the .star file as the group number for
% each particle.  Combine multiple successive micrographs if the group is
% too small.

doLoadFiles=1;  % option for skipping the long process.

if doLoadFiles
    disp('Getting a data.star file from the classification.');
    [stName, stPath]=uigetfile('*data.star','Select data star file');
    if isnumeric(stPath)  % user has clicked Cancel
        return
    end;
%     Pick up the local path string for display
    [basePath,starDir]=ParsePath(stPath);
    cd(basePath);  % cd to the enclosing folder
    [ba2,pa2]=ParsePath(basePath);
    [ba3,pa3]=ParsePath(ba2);
    localPath=[pa3 pa2 starDir];
    %%
    disp(['Reading ' stName '...']);
    [stBlocks,stData,ok]=ReadStarFile([starDir stName]);
    disp('done.');
    %%
    dat=stData{1};
    nim=numel(dat.rlnDefocusU);
    if isfield(dat,'rlnGroupNumber')
        minGroup=min(dat.rlnGroupNumber);
        maxGroup=max(dat.rlnGroupNumber);
        disp(['Min, max group values: ' num2str([minGroup maxGroup])]);
    else
        disp('No group field found.');
    end;

    ind=zeros(nim,1);
    for i=1:nim
%         ind(i)=rlDecodeImageName(dat.rlnImageName{i});
        ind(i)=sscanf(dat.rlnImageName{i},'%d');
    end;

    disp('Getting a si.mat file from the original images');
    [siName, siPath]=uigetfile('*si.mat','Select stack info file');
    if isnumeric(stPath)  % user has clicked Cancel
        return
    end;
    
    load([AddSlash(siPath) siName]);
end;
%% Get the micrograph number for each particle
% ------parameters---------
minNum=20;
doPickClasses=0;
desiredClasses=[3 4];

if doPickClasses
    clsOk=false(nim,1);
    for j=1:numel(desiredClasses)
        clsOk= (dat.rlnClassNumber==desiredClasses(j)) | clsOk;
    end;
    disp(['Class selection ' num2str(sum(clsOk)) ' out of ' num2str(numel(clsOk))]);
else
    clsOk=true(nim,1);
    disp('Using all classes');
end;
goodInds=ind(clsOk);
goodRows=find(clsOk);
ngim=sum(clsOk);  % number of good images

nmi=numel(si.mi);
def=zeros(nmi,1);
for i=1:nmi
    def(i)=si.mi{i}.ctf(1).defocus;  % defocus in um
end;

ddef=diff(def);

gNum=si.miIndex(goodInds); % starting group numbers


%% Figure out composite groups
h=hist(gNum,1:nmi);  % how many particles in each micrograph

% For each group containing at least one particle, search neighbors to find
% enough to add up to minNum
k=0;
group={[]};
i=1;
while i<=nmi  % index of micrographs
   if h(i)>0
       k=k+1;
       group{k}=i; % initialize a group
       i=i+1;
%        Grow the group
       while sum(h(group{k}))<minNum && i<nmi
           if h(i)>0
               group{k}=[group{k} i];
           end;
           i=i+1;
       end;  % group is done.
    else
       i=i+1;  % skip empty micrographs
   end;
end;
disp([num2str(k) ' groups defined from ' num2str(nmi)]);

% Assign group numbers
newGroupNumbers=zeros(ngim,1);
for i=1:k
    for j=1:numel(group{i})
        newGroupNumbers(gNum==group{i}(j))=i;
    end;
end;

% Insert the new group numbers into the old dat structure
dat.rlnGroupNumber=zeros(nim,1,'int16');
dat.rlnGroupNumber(goodRows)=newGroupNumbers;

% Let's sort by particle number
[sortInds,iMap]=sort(goodInds);

% Create the new structure containing only good rows, sorted by particle
newDat=struct;
fNames=fieldnames(dat);
for i=1:numel(fNames)
    field=fNames{i};
    q=dat.(field)(goodRows);
    newDat.(field)=q(iMap);
end;

% Save the new star file
[pa,nm,ex]=fileparts(stName); % construct the new name
p=strfind(nm,'data');
if numel(p)<1
    p=1;
end;
if doPickClasses
    str=sprintf('%d_', desiredClasses);
else
    str='ort_';
end;
outName=[starDir nm(1:p-1) 'cls' str nm(p:end) ex];
disp(['Writing ' outName]);  % write it
WriteStarFileStruct(newDat,'',outName);
disp('done.');


