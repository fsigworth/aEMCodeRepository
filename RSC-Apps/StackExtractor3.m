% StackExtractor3.m
% Quick, simplified version that works directly from the merged image.
% No padding or rotation; no CTF computation, no merge compensation.
% Stores *tstack.mrc
% and *tustack.mrc, the merged subtracted and unsubtracted stacks from all
% images, and *tsi.mat which contains the total si structure. Contrast is
% reversed so protein is white.
% Supports mi.particle.picks(:,10) flag for active particles.

batchMode=1;
allMiName='Info/allMis1.mat';
doDisplay=0;

boxSize=128;  % Size of boxes to be extracted from merged images.
%  boxSize=192;  % Size of boxes to be extracted from merged images.
% boxSize=256;
ds=2;        % downsampling of boxed particles from original micrograph
dfc=.1;      % Gauss filter for display (relative to original micrograph)

types=[16 32]; % flags for valid particles
restoreFromSiFile=0; % Use info from an si file instead of mi files.
loadAllMisFile=1;

% upsample the vesicle models?
upsampleSubtractedImage=0;
%altBasePath='../KvLipo121_2retrack/';
altBasePath='';
vus=2;  % default upsampling factor (used if we construct vesicle models)

minImageSD=.01;  % if std of image is below this, don't use it

doExtractParticles=1;
showAllParticles=0;  % display each particle as the program runs.
cpe0=0.8;     % default value, if not already in mi file
resetBasePath=0;
updateMiFile=0; % write out individual, modified mi files. Use this if you
% want individual mi files written if you use loadAllMisFile option, as
% well.

writeUnsubtractedStack=1;

tempStackName='TempStack2.mrc';
tempUStackName='TempUStack2.mrc';

% setOk4=0;  % force the unused field mi.vesicle.ok(:,4) to 1.
% mergeSuffixes={'' 'sf' 'su' 'si'}; % Suffix for merged image name, indexed by mergeMode.
inputModeSuffix=''; % expected suffix for input files

% Output file naming
stackDir='Stack2/';
dirVesicles='Vesicles/'; % location of modeled vesicle images
stackSuffix='tstack.mrc';
ustackSuffix='tustack.mrc';

% Output name is constructed thusly:
% [baseName weightString sizeString modeSuffix];
% e.g. name1w11p64m3
weightString=''; % string in output file name

defaultMembraneOffsetA=52;

% Where do we get vesicle models?
readVesicleModels=0;  % load vesicle images from files rather than making models on the fly.
readSubtractedImages=1; % load subtracted images from *mv files in Merged/ directory

writeVesicleModels=0;
writeSubtractedImages=0;
ignoreOldActiveFlags=0;

clear log % We had the error that the log function was overloaded.


% % %%%%%%%%%% Special options for non-particle stack %%%%%%%%
% boxSize=80;
% fHP=0;
% usePWFilter=0;
% types=[48 48];
% mergeMode=13;
% restoreFromSiFile=0;
% noDamage=0;
% setWeights=0;
% weightString='blanks';
% writeUnsubtractedStack=0;
% writeSubtractedImages=0;
% readSubtractedImages=0;
% % %%%%%%%%%%%%%  End special options


% dsScale=ds^2/4;  % fudge for vesicle model, don't really understand yet.
dsScale=1;   % scaling that depends on downsampling ratio
dds=2;       % further downsampling for micrograph display
vindex=0;    % 0 means force all vesicles to be modeled

if batchMode
    load(allMiName);
    nmi=numel(allMis);
    pa=[AddSlash(pwd) 'Info/'];
elseif restoreFromSiFile
    disp('Select si file');
    [oldSiName,pa]=uigetfile('*si.mat','Select si file to read');
    if isnumeric(pa)
        return
    end;
    oldSi=load([AddSlash(pa) oldSiName]);
    allMis=oldSi.si.mi;
    nmi=numel(allMis);
elseif loadAllMisFile
    disp('Select allMis.mat file');
    [allMisName,pa]=uigetfile('*.mat','Select allMis file to read');
    if isnumeric(pa)
        return
    end;
    load([AddSlash(pa) allMisName]);
    nmi=numel(allMis);
else % load individual mi files
    oldSiName='';
    disp('Select mi files');
    [miNames, pa]=uigetfile('*mi.txt','Select mi files','multiselect','on');
    if isnumeric(pa)
        return
    end;
    if ~iscell(miNames)
        miNames={miNames};
    end;
    nmi=numel(miNames);
end;

[rootPath, infoPath]=ParsePath(pa); % back out of Info or Stack directory
cd(rootPath);
disp(rootPath);

if ~exist(stackDir,'dir')
    mkdir(stackDir);
end;

%%
% Initialize the si structure
totalNParts=0;
np=1e6;  % preliminary stack size
si=struct;
si.miIndex=     uint16(zeros(np,1));
si.miParticle=  uint16(zeros(np,1));
si.alpha0=      single(zeros(np,1));
si.yClick=      single(zeros(np,1));  % in units of si.pixA
si.rVesicle=    single(zeros(np,1));
si.sVesicle=    single(zeros(np,1));
si.mi=cell(0);
nfiles=nmi;
badFlags=false(np,1);
iCoordsList=zeros(np,2,'single');
% ctfs=zeros(boxSize,boxSize,nfiles,'single');
% ctf1s=zeros(boxSize,boxSize,nfiles,'single');
% pwfs=zeros(boxSize,boxSize,nfiles,'single');
pixA0=0;  % unassigned value
fh=0;   % temp stack file handles
fhU=0;  % temp unsubtracted stack file
sumTotalStack=0;
sumTotalStackU=0;
if doDisplay
    figure(1);
    clf;
end;

fileIndex=1;
%%  Scan over files
while fileIndex<= nmi
    if restoreFromSiFile || loadAllMisFile
        mi=allMis{fileIndex};
%         disp(mi.baseFilename);
    else
        disp(['Reading ' miNames{fileIndex}]);
        mi=ReadMiFile([infoPath miNames{fileIndex}]);  % Load the mi file
    end;
    if ~isa(mi,'struct')  % a valid structure
        disp([num2str(fileIndex) ' mi is not a proper structure']);
        fileIndex=fileIndex+1;
        continue;
    end;
    mi.basePath=rootPath;  % Reassign this
    % Make new entries into the mi file
    mi.boxSize=boxSize;
    mi.stackPath=AddSlash(stackDir);
    si.mi{fileIndex}=mi;  % store a copy of the micrograph info even if there are
%                           no particles.
    
    if isfield(mi.particle,'picks') && numel(mi.particle.picks)>0
        if size(mi.particle.picks,2)<10 % don't have the flag field
            if fileIndex==1
                disp('Setting default pick(10) flags.');
            end;
            mi.particle.picks(:,10)=(typeArray>=types(1) & typeArray<=types(2)); % set them all active
            si.mi{fileIndex}=mi;  % update the copy of micrograph info
        end;
        active=mi.particle.picks(:,10)>0;
        nParts=sum(active);
        ourPicks=mi.particle.picks(active,:);
        ourMiParticles=find(active);
        typeArray=ourPicks(:,3);
    else
        nParts=0;
        ourPicks=[];
        disp([num2str(fileIndex) ' no picks.']);
        fileIndex=fileIndex+1;
        continue;
    end;
    n0=mi.imageSize/ds;
    %     Get the final pixel size
    pixA=mi.pixA*ds;
    if pixA0==0
        pixA0=pixA;
    end;
    if abs(pixA0-pixA)>.001
        warning(['Change in pixA values: ' num2str([pixA0 pixA]) '  ' miNames{fileIndex}]);
    end;
    if fileIndex==1
        disp(['   box size in A, box size in pixels: ' num2str([boxSize*pixA boxSize])]);
    end;
    
    disp(['Read images ' num2str(fileIndex) ' ' mi.baseFilename 'mi.txt  ' num2str(nParts)]);
    %     [mMergeU,pa,mImageOk]=meReadMergedImage(mi,0,inputModeSuffix);  % merged image
    %     subplot(1,2,1);
    %     imags(BinImage(mMergeU,dds));
    %     title(mi.baseFilename,'interpreter','none');
    %     drawnow;
    
    %         Read images, subtract vesicles
    %     if mImageOk && (nParts>0 || writeVesicleModels) % there are particles
    %             Read the images
    %       Scale the image by sqrt(doses(1)) to make it nominally unit variance.
    %       Also scale up by pixel size.
%     normScale=sqrt(mi.doses(1))*pixA/ds^2;  % not sure where factor of ds^2 comes from.
%             Read a vesicle image *v.mrc or a subtracted image *mv.mrc
%     vName=[dirVesicles mi.baseFilename 'v' inputModeSuffix '.mrc'];
%     if upsampleSubtractedImage
%         pa=[mi.basePath altBasePath];
%     else
%         pa=mi.basePath;
%         vus=1;  % no vesicle upsampling
%     end;
    mvName=[mi.procPath mi.baseFilename 'mv' inputModeSuffix '.mrc'];
    mvImageOk=exist(mvName,'file');
    mName=[mi.procPath mi.baseFilename 'm' inputModeSuffix '.mrc'];
    mImageOk=exist(mName,'file');
    imagesOk=mImageOk & mvImageOk;
    if imagesOk
%         disp(mvName);
        mvMergeU=ReadEMFile(mvName);
        mMergeU=ReadEMFile(mName);
        n0=size(mMergeU);
        n1=n0/ds;
        if any(size(mvMergeU)~=size(mMergeU))
            disp('image size discrepancy.');
            imagesOk=0;
        else
            mMerge=Downsample(mMergeU,n1)*ds^2;
            mvMerge=Downsample(mvMergeU,n1)*ds^2;
        end;
    end;
    imagesOk=imagesOk && std(mvMerge(:))>minImageSD && std(mMerge(:))>minImageSD;
    if ~imagesOk
        disp([num2str(fileIndex) ' missing or bad image.']);
        fileIndex=fileIndex+1;
        continue;
    end;
    
    %             We now have unfiltered mMergeU and mvMergeU.
    if doDisplay
        subplot(121);
    imags(GaussFilt(mvMerge,dfc*ds));
    title([num2str(fileIndex) ' ' mi.baseFilename],'interpreter','none');
    drawnow;
    end;
   
    %             Extract from a single image pair
    if doExtractParticles
        
        %                 msk=~meGetMask(mi,n0);  % Get the image mask
        %                 expandedMsk=msk;
        
        %
%         nPickerEntries=nParts;
        %                 if totalNParts==0  % the first time, allocate part of the total stack.
        %                     totalStack=single(zeros(boxSize,boxSize,nPicks));
        %                     totalStackU=single(zeros(boxSize,boxSize,nPicks));
        %                 end;
        nv=numel(mi.vesicle.x);
        vCoords=[mi.vesicle.x(:) mi.vesicle.y(:)];
        vRadii=mi.vesicle.r(:,1);

        stackSub=single(zeros(boxSize,boxSize,nParts));
        stackImg=single(zeros(boxSize,boxSize,nParts));
        sumImg=zeros(boxSize,boxSize);
        sumSub=zeros(boxSize,boxSize);
%         disp('Extracting');        
        %                 Particle extraction loop
        for i=1:nParts  % scan the particle coordinates
            type=ourPicks(i,3);
            vInd=ourPicks(i,4);
            
            pCoords=ourPicks(i,1:2);
            
            if vInd<1 && type<48  % no vesicle marked.  Find the nearest vesicle
                dists=sqrt(sum((vCoords-repmat(pCoords,nv,1)).^2,2));
                [minDist, vInd]=min(dists);
                if minDist > vRadii(vInd) % outside of membrane, check distance from mbn.
                    distm=abs(dists-vRadii);
                    [minDist, vInd]=min(distm);
                end;
                ourPicks(i,4)=vInd;  % insert our estimated vesicle index.
                disp(['Inserted vesicle index for particle, vesicle ' num2str([i vInd])]);
            end;
            %   Up to this point, all coordinates are in original pixel size.
            %   Make downsampled coordinates and do extraction
            iCoords=round(pCoords/ds)+1;  % downsampled coordinates
            %                 Store for box display
            iCoordsList(i,:)=iCoords;  % store for box display
            intCoords=round(iCoords);  % downsampled coordinates
%             fraCoords=iCoords-intCoords;
            xm=ExtractImage(mMerge,intCoords,boxSize);
            xmv=ExtractImage(mvMerge,intCoords,boxSize);
            %    Pad and reverse the contrast.
            ximg2=-xm;
            xsub2=-xmv;
            %                        Compute the angle relative to the y axis
            %                           positive alpha means particle lies ccw from vertical
            if vInd>0               % associated with a vesicle?
                pVec=pCoords-vCoords(vInd,:);
                alpha=atan2(pVec(2),pVec(1))-pi/2;
            else
                pVec=[0 0];
                alpha=0;
            end;
            stackSub(:,:,i)=xsub2;
            sumSub=sumSub+xsub2;
            
            stackImg(:,:,i)=ximg2;
            sumImg=sumImg+ximg2;
            
            if showAllParticles || (doDisplay && i==1)
                subplot(2,4,3);
                imags(xsub2);
                title(mi.baseFilename,'interpreter','none');
                %                         Show the particle image, and subtracted image
                subplot(2,4,4);
                imags(ximg2);
            end;
            
            if doDisplay && showAllParticles
                %                         Show the sums fsrom this micrograph
                subplot(2,4,7)
                imags(sumSub);
                title(i);
                subplot(2,4,8);
                imags(sumImg);
                drawnow;
            end;
            
            jt=i+totalNParts;  % pointer to overall stack index
            si.miIndex(jt)=     fileIndex;
            si.miParticle(jt)=  ourMiParticles(i); %% bug fixed here!!!
            si.alpha0(jt)=      alpha*180/pi;
            si.yClick(jt)=      hypot(pVec(1),pVec(2))/ds;  % in units of si.pixA
            if vInd>0
                si.rVesicle(jt)=    mi.vesicle.r(vInd)/ds;
                si.sVesicle(jt)=    mi.vesicle.s(vInd);
            end;
        end; % for i
        
        if fh==0 && nParts>0 % Set up to write out
            disp('Opening output files.');
            fh=WriteMRCHeader(stackImg,pixA,tempStackName);  % file handles for stack output files
            fhU=WriteMRCHeader(stackImg,pixA,tempUStackName);
        end;
        %                 Write the stack data
        fwrite(fh,stackSub,'float32');
        fwrite(fhU,stackImg,'float32');
        
        sumTotalStack=sumTotalStack+sum(stackSub,3);
        sumTotalStackU=sumTotalStackU+sum(stackImg,3);
        
        totalNParts=totalNParts+nParts;
        %                 figure(1);
        if doDisplay
        subplot(2,4,7)
        imags(sumTotalStack);  % unrotated image
        title(totalNParts);
        subplot(2,4,8);
        imags(sumTotalStackU);
        drawnow;
        end;
        
        if ~restoreFromSiFile && updateMiFile
            WriteMiFile(mi,[infoPath miNames{fileIndex}]);
            disp([infoPath miNames{fileIndex} ' updated']);
        end;
        
        %                 if writeIndividualStacks
        %                     outname=[mi.stackPath mi.baseFilename 'st.mrc'];
        %                     WriteMRC(stackSub,pixA,outname);
        %                     outname=[mi.stackPath mi.baseFilename 'stu.mrc'];
        %                     WriteMRC(stackImg,pixA,outname);
        %                 end;
    end; % if doExtractParticles
    fileIndex=fileIndex+1;
end; % while fileIndex
%
if fh>0 % we wrote something, finalize the files.
    frewind(fh);
    WriteMRCHeader(xsub2,pixA,tempStackName,[boxSize boxSize totalNParts],0,2,0,fh);
    fclose(fh);
    frewind(fhU);
    WriteMRCHeader(ximg2,pixA,tempUStackName,[boxSize boxSize totalNParts],0,2,0,fhU);
    fclose(fhU);
end;


if totalNParts>0 || ~doExtractParticles % We'll write out .si and perhaps stack data.
    %     update the si structure
    si.pixA=pixA;
    si.fHighPass=0;
    mi1=si.mi{1};
    if isfield(mi1,'ppVars')
        si.mbnOffset=mi1.ppVars.membraneOffsetA/pixA;
    else
        si.mbnOffset=defaultMembraneOffsetA/pixA;
    end;
    si.weights=mi.weights;
    np=totalNParts;
    % truncate the si arrays.
    si.miIndex=     si.miIndex(1:np,1);
    si.miParticle=  si.miParticle(1:np,1);
    si.alpha0=      si.alpha0(1:np,1);
    si.yClick=      si.yClick(1:np,1);
    si.rVesicle=    si.rVesicle(1:np,1);
    si.sVesicle=    si.sVesicle(1:np,1);
    si.ctfs = [];
    si.ctf1s=[];
    si.pwfs=[];
    if restoreFromSiFile && exist('oldActiveFlags','var')
        si.activeFlags=oldActiveFlags;
        si.activeFlagLog=oldActiveFlagLog;
    else
        si.activeFlags=true(np,1);
        si.activeFlagLog={[date '  StackExtractor3']};
    end;
    si.mergeMode=mi.mergeMode;
    %
    %%  Construct the stack file name
    %     Basename is up to 2nd underscore
    p=strfind(mi.baseFilename,'_');
    if numel(p)>1
        baseName=mi.baseFilename(1:p(2)-1);
    else
        baseName=mi.baseFilename;
    end;
    
    sizeString=sprintf('p%d',boxSize);
    modeSuffix='';
    if mi.mergeMode>1
        modeSuffix=['m' num2str(mergeMode)];
    end;
    
    outname=[mi.stackPath baseName sizeString modeSuffix];
    disp(['Base output name: ' outname]);
    siName=[outname 'tsi.mat'];
    if exist(siName,'file') && ~batchMode
        ch=MyInput(['Overwrite the file ' siName],'n');
        if ch~='y'
            outname=MyInput('New base name',[outname 'z']);
            siName=[outname 'tsi.mat']
        end;
    end;
    
    % save the stackInfo structure
    save([outname 'tsi.mat'],'si','-v7.3');
    disp(['Wrote the total stack info ' outname '.tsi.mat']);
    
    if totalNParts>0 % Actually particles to store, write the data.
        stackName=[outname stackSuffix];
        [status,result]=system(['mv ' tempStackName ' ' stackName]);
        disp(result);
        disp(['Wrote the stack ' outname stackSuffix]);
        
        ustackName=[outname ustackSuffix];
        [status,result]=system(['mv ' tempUStackName ' ' ustackName]);
        disp(result);
        disp(['Wrote the stack ' outname ustackSuffix]);
    else
        disp('No particle stacks written');
    end;
else
    disp('--nothing written.');
    disp(' ');
end;
