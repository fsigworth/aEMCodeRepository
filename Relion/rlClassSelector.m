% rlClassSelector
% From a data.star file from a classification step, pick one class and
% write out the new star file.  Also allows selection of Iso/Rso particles.

classSel=3;  % can be a vector of allowed classes
selectSided=1;
selectRso=1;

doLoadFiles=1;

if ~exist('doLoadFiles','var') || doLoadFiles
    disp('Getting a data.star file from the classification.');
    [stName, stPath]=uigetfile('*.star','Select data star file');
    if isnumeric(stPath)  % user has clicked Cancel
        return
    else
        cd(stPath);
    end
    [ba1,pa1]=ParsePath(stPath);
    [ba2,pa2]=ParsePath(ba1);
    [ba3,pa3]=ParsePath(ba2);
    localPath=[pa3 pa2 pa1];
    %%
    disp(['Reading ' stName '...']);
    [stBlocks,stData,ok]=ReadStarFile(stName);
    disp('done.');
    %%
    dat=stData{1};
    nim=numel(dat.rlnDefocusU);
    if isfield(dat,'rlnClassNumber')
        cls=dat.rlnClassNumber;
    else
        cls=zeros(nim,1);
    end
    def=(dat.rlnDefocusU + dat.rlnDefocusV)/2e4;  % defocus in um
    pInds=zeros(nim,1);
    disp('Decoding image numbers...');
    for i=1:nim
        pInds(i)=rlDecodeImageName(dat.rlnImageName{i});
    end
    disp('done.');
    
    select= any(dat.rlnClassNumber==classSel);
    disp([num2str(sum(select)) ' particles selected from class ' num2str(classSel)]);    
    
    if selectSided
        disp('Getting a si.mat file from the original images');
        [siName, siPath]=uigetfile('*si.mat','Select stack info file');
        if isnumeric(stPath)  % user has clicked Cancel
            return
        end;
        
        load([AddSlash(siPath) siName]);
        nsp=size(si.miIndex,1);  % number of si particles
        
        rso=zeros(nsp,1,'single');
        
        for i=1:nsp
            rso(i)=si.mi{si.miIndex(i)}.particle.picks(si.miParticle(i),7);
        end;
        
        select = select & (rso(pInds)==selectRso);
        disp([num2str(sum(select)) ' particles selected with rso=' num2str(selectRso)]);
    end;
%% 

disp('Getting the output star file name');
 [outName, outPath]=uiputfile('*.star');
 
%  sdat=reSplitStructureFields(dat,select);
 
 WriteStarFileStruct(dat,'',[outPath outName],select);
%% 
end;
