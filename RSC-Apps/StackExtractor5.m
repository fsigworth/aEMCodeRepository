% function StackExtractor5(names,pars)
nargin=0;
% Derived from StackExtractor4, but handles rawMicrographScaling and looks
% in the micrograph folder.
% Quick, simplified version that works directly from the merged image.
% No padding or rotation; no CTF computation, no merge compensation.
% Stores *tstack.mrc
% and *tustack.mrc, the merged subtracted and unsubtracted stacks from all
% images, and *tsi.mat which contains the total si structure. Contrast is
% reversed so protein is white.
% Supports mi.particle.picks(:,10) flag for active particles.
if nargin<2
    pars=struct; % Accept all the defaults
end;

dpars.batchMode= nargin>0 && numel(names)>0;

dpars.doDisplay=1;
dpars.doPrint=1; % print filename, number of particles to command window.
dpars.setAllActive=1;  % all valid particles are set to be active.
dpars.fHighPass=.002;  % A^-1. 0 means no highpass

dpars.boxSize=160;  % Size of boxes to be extracted from downsampled merged images.
dpars.ds=2;        % downsampling of images from which particles are boxed.
dpars.infoPath='Info/';
dpars.baseName='st';
dpars.outputDir='Stack1/';
dpars.inputModeSuffix=''; % set to 's' to use downsampled images.

dfc=.1;      % Gauss filter for display (relative to original micrograph)
types=[16 32]; % flags for valid particles
% how to get filenames (default is file selector)
dpars.restoreFromSiFile=0; % Use info from an si file instead of mi files.
loadAllMisFile=0;

pars=SetDefaultValues(dpars,pars,1); % 1 means check the fieldnames.

boxSize=pars.boxSize;

% upsample the vesicle models?
upsampleSubtractedImage=0;
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

tempStackName='TempStack2.mrc';
tempUStackName='TempUStack2.mrc';
stackSuffix='tstack.mrc';
ustackSuffix='tustack.mrc';

% setOk4=0;  % force the unused field mi.vesicle.ok(:,4) to 1.
% mergeSuffixes={'' 'sf' 'su' 'si'}; % Suffix for merged image name, indexed by mergeMode.

% Output name is constructed thusly:
% [pars.baseName weightString sizeString modeSuffix];
% e.g. name1w11p64m3
weightString=''; % string in output file name

defaultMembraneOffsetA=52;

clear log % We had the error that the log function was overloaded.

if pars.batchMode
    miNames=f2FindInfoFiles(dpars.infoPath);
    nmi=numel(miNames);
    pa=pars.infoPath;
elseif pars.restoreFromSiFile
    disp('Select si file');
    [oldSiName,pa]=uigetfile('*si.mat','Select si file to read');
    if isnumeric(pa)
        return
    end;
    oldSi=load([AddSlash(pa) oldSiName]);
    allMis=oldSi.si.mi;
    nmi=numel(allMis);
else % load individual mi files using file selector
    disp('Select mi files');
    [miNames, pa]=uigetfile('*mi.txt','Select mi files','multiselect','on');
    if isnumeric(pa)
        return
    end;
    if ~iscell(miNames)
        miNames={miNames};
    end;
    nmi=numel(miNames);
    for i=1:nmi
        miNames{i}=[AddSlash(pa) miNames{i}];
    end;
end;

[rootPath, infoPath]=ParsePath(pa); % back out of Info or Stack directory
cd(rootPath);
disp(rootPath);

CheckAndMakeDir(pars.outputDir,1);


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
iCoordsList=zeros(np,2,'single');
% ctfs=zeros(boxSize,boxSize,nfiles,'single');
% ctf1s=zeros(boxSize,boxSize,nfiles,'single');
% pwfs=zeros(boxSize,boxSize,nfiles,'single');
pixA0=0;  % unassigned value
fh=0;   % temp stack file handles
fhU=0;  % temp unsubtracted stack file
sumTotalStack=0;
sumTotalStackU=0;
if pars.doDisplay
    figure(1);
    clf;
end;

fileIndex=1;
%%  Scan over files
while fileIndex<= nmi
    if pars.restoreFromSiFile || loadAllMisFile
        mi=allMis{fileIndex};
    else
        if pars.doPrint
            disp(miNames{fileIndex});
        end;
        mi=ReadMiFile([miNames{fileIndex}]);  % Load the mi file
    end;
    if ~isa(mi,'struct')  % a valid structure
        disp([num2str(fileIndex) ' mi is not a proper structure']);
        fileIndex=fileIndex+1;
        continue;
    end;
    mi.basePath=rootPath;  % Reassign this
    
    % Make new entries into the mi file
    mi.boxSize=boxSize;
    mi.stackPath=AddSlash(pars.outputDir);
    si.mi{fileIndex}=mi;  % store a copy of the micrograph info even if there are
    %                           no particles.
    
    if isfield(mi.particle,'picks') && numel(mi.particle.picks)>0 % There are particles
        if size(mi.particle.picks,2)<10 || pars.setAllActive % don't have the flag field
            flags=mi.particle.picks(:,3);
            mi.particle.picks(:,10)=(flags>=types(1)) & (flags <=types(2)); % all valid particles are active
            si.mi{fileIndex}=mi;  % update the copy of micrograph info
        end;
        active=mi.particle.picks(:,10)>0;
        nParts=sum(active);
        ourPicks=mi.particle.picks(active,:);
        ourMiParticles=find(active);
        if pars.doPrint
            disp(nParts);
        end;
    else
        if pars.doPrint
            disp(['fileIndex ' num2str(fileIndex) ': no picks.']);
        end;
        fileIndex=fileIndex+1;
        continue;
    end;
    n0=mi.imageSize/pars.ds;
    %     Get the final pixel size
    pixA=mi.pixA*pars.ds;
    if pixA0==0
        pixA0=pixA;
    end;
    if abs(pixA0-pixA)>.01
        warning(['Change in pixA values: ' num2str([pixA0 pixA]) '  ' miNames{fileIndex}]);
    end;
    if fileIndex==1
        disp(['   box size in A, box size in pixels: ' num2str([boxSize*pixA boxSize])]);
    end;

    [mMergeU,mMat,mImageOk]=meLoadNormalizedImage(mi,mi.padImageSize);
    [mvMergeU,mvMat,mvImageOk]=meLoadNormalizedImage(mi,mi.padImageSize,'mv');
    
%     mvName=[mi.procPath mi.baseFilename 'mv' pars.inputModeSuffix '.mrc'];
%     mvImageOk=exist(mvName,'file');
%     mName=[mi.procPath mi.baseFilename 'm' pars.inputModeSuffix '.mrc'];
%     mImageOk=exist(mName,'file');
    imagesOk=mImageOk & mvImageOk;
    if imagesOk
        if isfield(mi,'useMicrographCoords') && mi.useMicrographCoords
            pOffset=-mMat(1:2,3)'; % particle padding offset
        else
            pOffset=[0 0];
        end
%         mvMergeU=ReadEMFile(mvName);
%         mMergeU=ReadEMFile(mName);
        n0=size(mMergeU);
        n1=n0/pars.ds;
        if any(size(mvMergeU)~=size(mMergeU))
            disp('image size discrepancy.');
            imagesOk=0;
        else
            mMerge=Downsample(mMergeU,n1)*pars.ds^2;
            mvMerge=Downsample(mvMergeU,n1)*pars.ds^2;
        end;
    end;
    imagesOk=imagesOk && std(mvMerge(:))>minImageSD && std(mMerge(:))>minImageSD;
    if ~imagesOk
        disp([num2str(fileIndex) ' missing or bad image.']);
        fileIndex=fileIndex+1;
        continue;
    end;
    
    %             We now have unfiltered mMerge and mvMerge, with pixA
    %             being the pixel size.
    if pars.fHighPass>0
        mvfMerge=GaussHP(mvMerge,pars.fHighPass*pixA);
        mfMerge=GaussHP(mMerge,pars.fHighPass*pixA);
    else
        mvfMerge=mvMerge;
        mfMerge=mMerge;
    end;
    if pars.doDisplay % Show the micrograph
        subplot(121);
        imags(GaussFilt(mvfMerge,dfc*pars.ds));
        title([num2str(fileIndex) '  ' mi.baseFilename '  ' num2str(nParts)],'interpreter','none');
        drawnow;
    end;
    
    % -----------Accumulate particles from a micrograph pair -------------
    if doExtractParticles
        nv=numel(mi.vesicle.x);
        vCoords=[mi.vesicle.x(:) mi.vesicle.y(:)];
        vRadii=mi.vesicle.r(:,1);
        
        stackSub=single(zeros(boxSize,boxSize,nParts));
        stackImg=single(zeros(boxSize,boxSize,nParts));
        sumImg=zeros(boxSize,boxSize);
        sumSub=zeros(boxSize,boxSize);
        %                 Particle extraction loop
        for i=1:nParts  % scan the particle coordinates
            type=ourPicks(i,3);
            vInd=ourPicks(i,4);
            
            pCoords=ourPicks(i,1:2)+pOffset;
            
            if vInd<1 && type<48  % no vesicle marked.  Find the nearest vesicle
                dists=sqrt(sum((vCoords-repmat(pCoords,nv,1)).^2,2));
                [minDist, vInd]=min(dists);
                if minDist > vRadii(vInd) % outside of membrane, check distance from mbn.
                    distm=abs(dists-vRadii);
                    [~, vInd]=min(distm);
                end;
                ourPicks(i,4)=vInd;  % insert our estimated vesicle index.
                disp(['Inserted vesicle index for particle, vesicle ' num2str([i vInd])]);
            end;
            %   Up to this point, all coordinates are in original pixel size.
            %   Make downsampled coordinates and do extraction
            iCoords=round(pCoords/pars.ds)+1;  % downsampled coordinates
            %                 Store for box display
            iCoordsList(i,:)=iCoords;  % store for box display
            intCoords=round(iCoords);  % downsampled coordinates
            xm=ExtractImage(mfMerge,intCoords,boxSize);
            xmv=ExtractImage(mvfMerge,intCoords,boxSize);
            
            %    reverse the contrast.
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
            
            stackSub(:,:,i)=xsub2; % insert the image into the stack
            sumSub=sumSub+xsub2;
            
            stackImg(:,:,i)=ximg2; % insert the unsubtracted image
            sumImg=sumImg+ximg2;
            
            if showAllParticles || (pars.doDisplay && i==1) % show the first particle
                % in the micrograph
                subplot(2,4,3);
                imags(xsub2);
                %                 title(mi.baseFilename,'interpreter','none');
                %                         Show the particle image, and subtracted image
                subplot(2,4,4);
                imags(ximg2);
            end;
            
            if pars.doDisplay && showAllParticles
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
            si.yClick(jt)=      hypot(pVec(1),pVec(2))/pars.ds;  % in units of si.pixA
            if vInd>0
                si.rVesicle(jt)=    mi.vesicle.r(vInd)/pars.ds;
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
        if pars.doDisplay
            subplot(2,4,7)
            imags(sumTotalStack);  % unrotated image
            title(totalNParts);
            subplot(2,4,8);
            imags(sumTotalStackU);
            drawnow;
        end;
        
        if ~pars.restoreFromSiFile && updateMiFile
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
disp(['Total particles extracted: ' num2str(totalNParts)]);

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
    si.pars.fHighPass=0;
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
    if pars.restoreFromSiFile && exist('oldActiveFlags','var')
        si.activeFlags=oldActiveFlags;
        si.activeFlagLog=oldActiveFlagLog;
    else
        si.activeFlags=true(np,1);
        si.activeFlagLog={[date '  StackExtractor3']};
    end;
    if isfield(mi,'mergeMode')
        si.mergeMode=mi.mergeMode;
    else
        si.mergeMode=0;
    end;
    %
    %%  Construct the stack file name
    
    sizeString=sprintf('p%g',boxSize);
    modeSuffix='';
    if si.mergeMode>1
        modeSuffix=['m' num2str(si.mergeMode)];
    end;
    if pars.ds>1
        dsSuffix=sprintf('ds%g',pars.ds);
    end;
    if pars.fHighPass>0
        filtSuffix=sprintf('fh%g',round(1000*pars.fHighPass));
    end;
    % -------------Name constructed here------------
    outname=[mi.stackPath pars.baseName sizeString dsSuffix modeSuffix filtSuffix];
    disp(['Base output name: ' outname]);
    
    siName=[outname 'tsi.mat'];
    if exist(siName,'file') && ~pars.batchMode
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
