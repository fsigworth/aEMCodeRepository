% rsStackExtractClass.m
% After 3D clustering, select particles that belong to a particular class.

writeParticles=0;  % write a new stack from the class


disp('Select the roi.mat file.');
[roiName,roiPath]=uigetfile('*roi.mat','Select an roi.mat file');
if isnumeric(roiName)
    return
end;

load([roiPath roiName]);

disp('Loading the ri.mat file');
load([roiPath 'ri.mat']);

cd([roiPath '../../']);

%%
disp('Select the si file.');
[stackName,stackPath]=uigetfile('*si.mat','Select an si.mat file');
load([stackPath stackName]);
%%
nVols=size(roi.pVols,1);
nActive=size(roi.pVols,2);
disp([num2str(nActive) 'active particles']);
pVols=sum(roi.pVols,2)/nActive;
disp('Class  fraction');
for i=1:nVols
    disp([i pVols(i)]);
end;
iVol=MyInput('Volume to use',1);

% iFlags=roi.pVols(iVol,:)==max(roi.pVols);
% iFrac=MyInput('Percentile',99);
% iProbs=roi.pVols(iVol,iFlags);
% %
% val=Percentile(iProbs,iFrac/100)
%
pThresh=.7;
iFlags=roi.pVols(iVol,:)>pThresh;
disp(['Particles to save: ' num2str(sum(iFlags))]);

origPtrs=find(ri.activeFlags);
absPtrs=origPtrs(iFlags);
absFlags=false(size(ri.activeFlags));
absFlags(absPtrs)=true;

if writeParticles  % we'll store the truncated stack
        
    disp('Getting output file name.');
    [pa,nm,ex]=fileparts(stackName);
    baseStackName=nm(1:end-2);
    outName1=[baseStackName 'Class' num2str(iVol) nm(end-1:end) ex];
    [outStackName stackPath]=uiputfile(outName1);
    %%
    inputImagesName=[baseStackName 'stack.mrc'];
    disp(['Reading stack data ' inputImagesName '...']);
    [imgs,s]=ReadMRC([stackPath inputImagesName]);
    
    inputAltName=[baseStackName 'ustack.mrc'];
    if exist([stackPath inputAltName],'file')
        disp(['Reading stack data ' inputImagesName '...']);
        [altImgs,s]=ReadMRC([stackPath inputAltName]);
    end;
    disp(' Writing files.');
    %%
    % Change the variable names for saving
    si0=si;
    [si,imgss,altImgss]=rsStackSplit(absFlags,si0,imgs,altImgs);

    [pa,nm,ex]=fileparts(outStackName);
    baseOutputName=nm(1:end-2);

%     Save, then restore the old variable names
    save([stackPath baseOutputName 'si.mat'],'si');
    sis=si;
    si=si0;
    si0=[];

    WriteMRC(imgss,s.pixA,[stackPath baseOutputName 'stack.mrc']);
    WriteMRC(altImgss,s.pixA,[stackPath baseOutputName 'ustack.mrc']);
    disp('done.');

else
    nfe=size(si.activeFlags,2)+1; % number of flag entries
    si.activeFlags(:,nfe)=absFlags;
    [~,roiPath1]=ParsePath(roiPath);
    si.activeFlagLog{nfe,1}=[date '  rsStackExtractClass, class ' num2str(iVol)...
        ' using ' roiPath1 roiName];
    ok=input(['Overwrite ' stackName '? '],'s');
    if ok=='y'
        save([stackPath stackName],'si');
        disp([stackPath stackName ' written.']);
    else
        disp('nothing written.');
    end;
end;

