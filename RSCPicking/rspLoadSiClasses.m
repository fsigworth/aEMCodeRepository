% function [dis,si,miNames]=rspLoadSiClasses(dis);
% Load a star particles file from a relion Classification or Select job,
% and the corresponding si file. Update the si.mi.particle.picks to include
% class numbers. Return si and a list of mi filenames (needed but not actually used by
% SimpleRSPicker.)

% miNames={};
doSortClasses=1;


% % Pick up an si file.
    disp(['Getting an si file']);
     [siName,siPath]=uigetfile('*si.mat');
    if ischar(siPath)
        disp(['Loading ' siPath siName]);
        load([siPath siName]);
    else
        si=[];
        return
    end;
%
disp('Getting the corresponding star files');
[dataStarName,starPath]=uigetfile('*_data.star');

baseName=dataStarName(1:end-10);
modelStarName=[baseName '_model.star'];
disp(['Loading ' starPath modelStarName]);
[mdn,md]=ReadStarFile([starPath modelStarName],1);
pClasses=find(strcmp('data_model_classes',mdn),1);
md=md{pClasses}; % model class info
ncls=numel(md.rlnClassDistribution);
if doSortClasses
    [clsProbs,clsInds]=sort(md.rlnClassDistribution,'descend');
    clsInvInds=zeros(ncls,1);
    for i=1:nCls
        clsInvInds(i)=find(clsInds==i);
    end;
else
    clsProbs=md.rlnClassDistribution;
    clsInds=1:numel(clsProbs);
    clsInvInds=clsInds;
end;
%%
disp('Loading the class images');
mrcName=[baseName '_classes.mrcs'];
imgs=ReadMRC([starPath mrcName]);
sortImgs=imgs(:,:,clsInds);
ncls=size(imgs,3);
labels=cell(ncls,1);
for i=1:ncls
    labels{i}=sprintf('%4.2g%%    \\bf{%g}',100*clsProbs(i),i);
end;
figure(1);
imagsar(imgs,.001,2,labels);

disp(['Loading ' starPath dataStarName]);
[dn,d]=ReadStarFile([starPath dataStarName],1);
d=d{1};
nl=numel(d.rlnImageName);


%% Match the stack name to particle stack name to select the particle set.

disp('Decoding and matching particle references.');
siStackName=[siName(1:end-6) 'stack.mrcs']; % predicted stack file name
partStackFlag=false(nl,1);
partIndex=zeros(nl,1);
for i=1:nl
    str=d.rlnImageName{i};
    p=strfind(str,'@');
    if numel(p)==1
        partStackFlag(i)=strcmp(str(p+1:end),siStackName);
        partIndex(i)=str2double(str(1:p-1));
    end;
end;
npa=sum(partStackFlag);
disp([num2str(npa) ' particles matched.']);

if npa==0
%     si=[];
    return
end;

% Add the particle.class field to the relevant si.mi structures
disp('Inserting the mi.particle.class arrays.');
nmi=numel(si.mi);
for i=1:nmi
    if isfield(si.mi{i}.particle,'picks')
        np=size(si.mi{i}.particle.picks,1);
    else
        si.mi{i}.particle.picks=zeros(0,14,'single');
        np=0;
    end;
    si.mi{i}.particle.class=zeros(np,1,'single');
end;
activePartIndex=partIndex(partStackFlag);
activePartClass=d.rlnClassNumber(partStackFlag);
    miInd=si.miIndex(activePartIndex);
    paInd=si.miParticle(activePartIndex);
for i=1:npa
    si.mi{miInd(i)}.particle.class(paInd(i))=clsInvInds(activePartClass(i));
end;

% miNames=cell(nmi,1);
% for i=1:nmi
%     miNames{i}=[si.mi{i}.baseFilename 'mi.txt'];
% end;
% 
[~,nm]=fileparts(siName);
disp(['Saving ' AddSlash(pwd) nm '_classes.mat']);
save([nm '_classes.mat'],'si');

