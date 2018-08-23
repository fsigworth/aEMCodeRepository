% rlReplaceStackName
    origPath=pwd;
    disp('Getting the reference data.star file.');
    [stName, stPath]=uigetfile('*.star','Select data star file');
    if isnumeric(stPath)  % user has clicked Cancel
        return
    end;
    %%
    disp(['Reading ' stName '...']);
    [stBlocks,stData,ok]=ReadStarFile([AddSlash(stPath) stName]);
    %%
    disp('Decoding names.');
    dat=stData{1};
    nim=numel(dat.rlnImageName);
    ind=zeros(nim,1);
    files=cell(nim,1);
    for i=1:nim
       [ind(i),files{i}]=rlDecodeImageName(dat.rlnImageName{i});
    end;
%%
    disp(['First stack name is: ' files{1}]);
    
    disp('Getting the new data directory');
    datPath=uigetdir;
    if isnumeric(datPath)
        return
    end;
 datPath=GetRelativePath(datPath)
 %   datPath=AddSlash(datPath)
    disp('Constructing new filenames');
    imgNames=cell(nim,1);
    for i=1:nim
        [pa,nm,ex]=fileparts(files{i});
        newName=[datPath nm ex];
        imgNames{i}=sprintf('%05d@%s',ind(i),newName);
    end;
    
    dat2=dat;
    dat2.rlnImageName=imgNames;
    
    disp('Getting the place to store the new particles.star.');
    outPath=uigetdir;
    outPath=GetRelativePath(outPath)
    disp('Writing...');
    WriteStarFileStruct(dat2,'',[outPath 'particles.star']);
    disp('done.');
 