% rlPruneParticles
% Purge lines from a .star file according to the inactive array.
% 
% doLoadFiles=1;
% 
% if ~exist('doLoadFiles','var') || doLoadFiles
    disp('Getting a data.star file from the classification.');
    [stName, stPath]=uigetfile('*.star','Select data star file');
    if isnumeric(stPath)  % user has clicked Cancel
        return
    else
        cd(stPath);
    end
    
    disp('Getting the output file name');
    [outName,outPath]=uiputfile('*.star','Output star file name');
    if isnumeric(outPath)
        return
    end;

    disp('Getting a mapping file');
    [matName, matPath]=uigetfile('*.mat','Select the mapping file');
    load([AddSlash(matPath) matName]);
    inactive=find(pInd2<1);
    disp([num2str(numel(inactive)) ' inactive out of ' num2str(numel(pInd2))]);


    %%
    disp(['Reading ' stName '...']);
    [stBlocks,stData,ok]=ReadStarFile(stName);
    disp('done.');
    %%
    dat=stData{1};
    nim=numel(dat.rlnDefocusU);
    pInds=zeros(nim,1);
    active=false(nim,1);
    for i=1:nim
        ind=rlDecodeImageName(dat.rlnImageName{i});
        active(i)=~any(inactive==ind);
    end
    disp([num2str(sum(~active)) ' inactive out of the set ' num2str(nim)]);
    
    WriteStarFileStruct(dat,'images',[AddSlash(outPath) outName],active);
    