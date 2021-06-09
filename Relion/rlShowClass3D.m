% rlShowClass3D.m

displayAngle=45;
% displayAngle=90
zFraction=.6; % relative Z-height of slice
skipLoading=1;

if ~(skipLoading && exist('maps','var'))
[modelName,modelPath]=uigetfile('*model.star');
% cd(modelPath);
modelName=[modelPath modelName];
[pa,nm,ex]=fileparts(modelName);
if ~strcmp(ex,'.mrc')
    p=strfind(modelName,'_');
    if numel(p)<2
        p(2)=numel(nm);
    end;
    baseName=modelName(1:p(2));

    [blocks,dat,ok]=ReadStarFile(modelName);
    if ~ok
        error('No data in Star file.');
    end;

    q1=strcmp('data_model_classes',blocks);
        if ~any(q1)
            error('No model_classes block found'); 
        end;
        classesIndex=find(q1,1);
        q1=strcmp('data_model_general',blocks);
        if ~any(q1)
            error('No model_general block found');
        end;
        generalIndex=find(q1,1);

    n=dat{generalIndex}.rlnOriginalImageSize;
    pixA=dat{generalIndex}.rlnPixelSize;
    %%
    mapNames=dat{classesIndex}.rlnReferenceImage;
    mapProbs=dat{classesIndex}.rlnClassDistribution;
    if ischar(mapNames)
        nim=1;
        mapNames={mapNames};
    else
        nim=numel(mapNames);
    end;
else
    nim=1;
    mapNames={modelName};
    mapProbs=1;
end;
%%
maps=zeros(n,n,n,nim,'single');
%
opts=struct;
opts.rowLabels=cell(nim,1);
for i=1:nim
    probString=num2str(mapProbs(i));
    disp([mapNames{i} '  ' probString]);
    opts.rowLabels{i}=probString;
    maps(:,:,:,i)=ReadMRC(mapNames{i});
end;
end; % if ~skipLoading
%%
n=size(maps,1);
n2=ceil((n+1)/2);
nz=round(n*zFraction);
ShowSections(maps,[n2,n2,nz],displayAngle); % not using opts because ShowSections can't make labels.

return;

%%
toShow=[1 4 7]
toShow=[2 5]
for i=toShow;
    figure(i);
    q=Downsample(maps(:,:,:,i),192);
    imagsar(q(:,:,1:85),.0001,1);
end;
