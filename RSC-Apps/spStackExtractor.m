% spStackExtractor
% Derived from StackExtractor2m.m
% Simplified version that works directly from the merged image.
% For each particle, extract its region from the merged image; rotate the
% merged particle to the standard position, if desired
% (writeIndividualStacks=1) write *st.mrc and *stu.mrc files (stack and
% unsubtracted stack) into the Stack/ directory.  Also store *tstack.mrc
% and *tustack.mrc, the merged subtracted and unsubtracted stacks from all
% images, and *tsi.mat which contains the total si structure. Contrast is
% reversed so protein is white.

rsMode=0;
boxSize=96;  % Size of boxes to be extracted from merged images.
% boxSize=256;  % Size of boxes to be extracted from merged images.
fHP=.001;    % highpass Gauss filter in A^-1
ds=2;        % downsampling of boxed particles from original micrograph
usePWFilter=1;  % prewhitening filter
doExtractParticles=1;
cpe0=16;     % default value, if not already in mi file
resetBasePath=1;
restoreFromSiFile=0; % Use info from an si file instead of mi files.
updateMiFile=0;
useCorrection=0;  % Use local vesicle modeling to refine vesicle subtraction
useCircularVesicleModels=0;  % ignore higher-order terms in the vesicle model, for experimentation purposes
weightString=''; % string in output file name
modeString='v2';

writeIndividualStacks=0;  % Write stack files for each micrograph too.
writeUnsubtractedStack=0;
setOk4=1;  % force the unused field mi.vesicle.ok(:,4) to 1.
stackSuffix='tstack.mrc';
ustackSuffix='tustack.mrc';
defaultMembraneOffsetA=52;
readVesicleModels=0;  % load vesicle images from files rather than making models on the fly.
readSubtractedImages=1; % load subtracted images from files
writeVesicleModels=0;
writeSubtractedImages=1;
nZeros=1;
stackDir='Stack/';
dirVesicles='Vesicles2/'; % location of modeled vesicle images

% dsScale=ds^2/4;  % fudge for vesicle model, don't really understand yet.
dsScale=1;
dds=2;       % further downsampling for micrograph display
vindex=0;    % 0 means force all vesicles to be modeled
padSize=NextNiceNumber(boxSize*1.5);
% Make the upsampled pad mask, and the particle fourier mask
padMask=fuzzymask(padSize,2,padSize*.48,padSize*.04);
fmasko=ifftshift(fuzzymask(padSize,2,.45*padSize,.05*padSize));

if restoreFromSiFile
    [oldSiName,pa]=uigetfile('*si.mat','Select si file to read');
    oldSi=load([AddSlash(pa) oldSiName]);
    miFiles=oldSi.si.mi;
    oldActiveFlags=oldSi.si.activeFlags;
    oldActiveFlagLog=oldSi.si.activeFlagLog;
    %     clear si;
else
    oldSiName='';
    [miFiles, pa]=uigetfile('*mi.txt','Select mi files','multiselect','on');
    if ~iscell(miFiles)
        miFiles={miFiles};
    end;
end;

if isnumeric(pa)  % user clicked Cancel
    return
end;

[rootPath, infoPath]=ParsePath(pa);
cd(rootPath);
if ~exist(stackDir,'dir');
    mkdir(stackDir);
end;

%%
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
nfiles=numel(miFiles);
badFlags=false(np,1);
iCoordsList=zeros(np,2,'single');
ctfs=zeros(boxSize,boxSize,nfiles,'single');
ctf1s=zeros(boxSize,boxSize,nfiles,'single');
pwfs=zeros(boxSize,boxSize,nfiles,'single');
pixA0=0;  % unassigned value

figure(1);
clf;
SetGrayscale;

fileIndex=1;
%%
while fileIndex<= numel(miFiles)
    if restoreFromSiFile
        mi=miFiles{fileIndex};
    else
        disp(['Reading ' miFiles{fileIndex}]);
        mi=ReadMiFile([infoPath miFiles{fileIndex}]);  % Load the mi file
    end;
    if isa(mi,'struct')  % a valid structure
        disp(mi.baseFilename)
        mi.basePath=rootPath;  % Reassign this
        if fileIndex==1
            modelSpectrum=CCDModelSpectrum2D(mi.camera);  % handle DE-12 or CCD
        end;
        % Make new entries into the mi file
        mi.boxSize=boxSize;
        if ~isfield(mi,'stackPath')
            mi.stackPath='';
        end;
        mi.stackPath=AddSlash(stackDir);
        if ~isfield(mi,'cpe')
            mi.cpe=cpe0;
        end;
        
        if isfield(mi.particle,'picks')
            nPicks=size(mi.particle.picks,1);
        else
            nPicks=0;
        end;
        n0=mi.imageSize/ds;
        %     Get the final pixel size
        pixA=mi.pixA*ds;
        if pixA0==0
            pixA0=pixA;
        end;
        if abs(pixA0-pixA)>.001
            warning(['Change in pixA values: ' num2str([pixA0 pixA]) '  ' miFiles{fileIndex}]);
        end;
        if fileIndex==1
            disp(['   box size in A / pixels: ' num2str([round(boxSize*pixA) boxSize])]);
        end;
        
        si.mi{fileIndex}=mi;  % store a copy of the micrograph info
        if nPicks>0 || writeVesicleModels % there are particles, insert filter functions
            
            if usePWFilter
                %     particle-sized noise-whitening filter
                nZeros=1;
                Tic=meGetNoiseWhiteningFilter(mi,boxSize,ds,nZeros,fHP*pixA);  % for particle ctf
            else
                Tic=1;
            end;
            %             Set up CTF and other filters for this micrograph
            freqs=RadiusNorm(boxSize)/pixA;
            [pcoeffs, pctf, pdctfs]=meComputeMergeCoeffs2(freqs,mi);
            
            ctfs(:,:,fileIndex)=pctf.*Tic;
            %             Accumulate the other filters in case we have to undo the
            %             merging.  To change an image to mode 2 merging,
            %             do this:
            %               ctfEff=ctfs./pwfs;
            %               filt=abs(ctf1s)./ctfEff;  % no abs for mode 3
            %               img=real(ifftn(fftn(img).*ifftshift(filt)));
            %               ctfs=ctf1s;
            ctf1s(:,:,fileIndex)=pdctfs(:,:,1); % first ctf
            pwfs(:,:,fileIndex)=Tic;  % pre-whitening filter that was applied
            
            %             if doExtractParticles  % go ahead and process micrographs
            if usePWFilter
                [Ti,filterOk]=meGetNoiseWhiteningFilter(mi,n0,ds,nZeros,fHP*pixA);  % for micrograph
                if ~filterOk
                    disp('No specimen-noise model, using simple highpass filter.');
                end;
            else
                Ti=1;
            end;
            disp('Read images');
            [mMergeU,pa,ok]=meReadMergedImage(mi,0);  % merged image
            %       Scale it up by sqrt(doses(1)) to make it nominally unit variance.
            mMerge=real(ifftn(fftn(mMergeU).*ifftshift(Ti))).*sqrt(mi.doses(1)); % prewhitened
            
            subplot(1,2,1);
            imacs(BinImage(mMerge,dds));
            title(mi.baseFilename,'interpreter','none');
            drawnow;
            
            vName=[dirVesicles mi.baseFilename 'v.mrc'];
            mvName=[mi.basePath mi.procPath mi.baseFilename 'mv.mrc'];
            [mvName,mvOk]=CheckForImageOrZTiff(mvName);
            indBadVesicles=[];
            if readSubtractedImages && mvOk
                mvMergeU=ReadEMFile(mvName);
                mvMerge=real(ifftn(fftn(mvMergeU).*ifftshift(Ti)))*sqrt(mi.doses(1));
            elseif readVesicleModels && exist(vName,'file')
                vMergeU=ReadMRC(vName);
                mvMergeU=mMergeU-vMergeU;
                vMerge=real(ifftn(fftn(vMergeU).*ifftshift(Ti)))*sqrt(mi.doses(1));
                mvMerge=mMerge-vMerge;
                %                         okName=['Vesicles/' mi.baseFilename 'ok.mat'];
                %                         sok=load(okName);
                %                         indBadVesicles=find(~sok.vOk);
            else
                disp('Compute vesicle models');
                vesMi=mi;
                if useCircularVesicleModels  % truncate to the zeroth term in r and s
                    vesMi.vesicle.r=real(mi.vesicle.r(:,1));
                    vesMi.vesicle.s=real(mi.vesicle.s(:,1));
                end;
                if setOk4
                    vesMi.vesicle.ok(:,4)=true;  % Ignore the 4th column of this.
                end;
                doCTF=1;
                %         We fit only the vesicles that have particles
                %         vindex=find(hist(mi.particle.picks(:,4),1:nPicks)>0);
                %         vSet=meMakeModelVesicleSet(mi,fitM,vindex,doCTF,Ti,nZeros)*dsScale;
                %           make vesicles with amplitude fitting
                vMergeU=meMakeModelVesicles(vesMi,size(mMerge),vindex,doCTF,0);
                vMerge0=real(ifftn(fftn(vMergeU).*ifftshift(Ti)))*sqrt(mi.doses(1));  % prewhitened
                %                 vSet0=meMakeModelVesicleSet(vesMi,mMerge,vindex,doCTF,Ti,nZeros,1);
                %                 vMerge0=sum(vSet0,3);
                vsCorr=zeros(size(vMerge0),'single');
                if useCorrection
                    imacs(GaussFilt(mvMerge0,.2));
                    title('before correction');
                    drawnow;
                    %         Get the residuals, which will be subtracted from individual
                    %         exposures.  the residual computation is a linear operation.
                    vsCorr=meGetVesicleResiduals(mi,mvMerge0);
                    vMerge=vMerge0+vsCorr;
                else
                    vMerge=vMerge0;
                end;
                mvMergeU=mMergeU-vMergeU;
                mvMerge=mMerge-vMerge;
                %                 if writeVesModels
                %                     subplot(1,2,2);
                %                     imags(mMerge-vMerge);
                %                     drawnow;
                %                     vesName=['CircVesModels/' mi.baseFilename 'ves.mrc'];
                %                     disp(['Writing ' vesName]);
                %                     WriteMRC(vMerge*mi.doses(1),pixA,vesName);
                %                     vesMatName=['CircVesModels/' mi.baseFilename 'ves.mat'];
                %                     disp(['Writing ' vesMatName]);
                %                     ves=vMerge*mi.doses(1);
                %                     save(vesMatName,'ves');
                %                 end;
                if writeSubtractedImages
                    [pa,nm,ex]=fileparts(mvName);
                    altName=[AddSlash(pa) nm 'z.tif'];  % compressed name
                    pars.snrRatio=200;
                    WriteZTiff(mvMergeU,pixA,altName,pars);
                    disp([nm 'z.tif written.']);
                end;
            end;
            if doExtractParticles
                nPicks=nPicks * (ok>0);  %%%% don't pick anything if no vesicle image.
                %                 disp('Compute vesicle models');
                %                 if setOk4
                %                     mi.vesicle.ok(:,4)=true;  % Ignore the 4th column of this.
                %                 end;
                %                 doCTF=1;
                %                 %         We fit only the vesicles that have particles
                %                 %         vindex=find(hist(mi.particle.picks(:,4),1:nPicks)>0);
                %                 fitM=size(mMerge);
                %                 %         vSet=meMakeModelVesicleSet(mi,fitM,vindex,doCTF,Ti,nZeros)*dsScale;
                %                 %           make vesicles with amplitude fitting
                %                 vSet0=meMakeModelVesicleSet(mi,mMerge,vindex,doCTF,Ti,nZeros,mergeMode);
                %                 vsCorr=zeros(size(vSet0),'single');
                %                 mvSet0=mSet-vSet0;  % subtracted images
                %                 if useCorrection
                %                     imacs(GaussFilt(sum(mvSet0,3),.2));
                %                     title('before correction');
                %                     drawnow;
                %                     %         Get the residuals, which will be subtracted from individual
                %                     %         exposures.  the residual computation is a linear operation.
                %                     for i=1:nim
                %                         vsCorr(:,:,i)=meGetVesicleResiduals(mi,sum(mvSet0,3));
                %                     end;
                %                     vSet=vSet0+vsCorr;
                %                 else
                %                     vSet=vSet0;
                %                 end;
                %
                %         imacs(BinImage(mvMerge,dds));
                imacs(GaussFilt(mvMerge,.2));
                title(mi.baseFilename,'interpreter','none');
                drawnow;
                % mtf=TestMarkParticles(mi,mMerge);
                
                % compute the CTF, including the prewhitening filter effects,
                % for individual particles
                freqs=RadiusNorm(boxSize)/pixA;
                [pcoeffs, pctf, pdctfs]=meComputeMergeCoeffs2(freqs,mi);
                
                ctfs(:,:,fileIndex)=pctf.*Tic;
                %             Accumulate the other filters in case we have to undo the
                %             merging
                ctf1s(:,:,fileIndex)=pdctfs(:,:,1); % first ctf
                pwfs(:,:,fileIndex)=Tic;  % pre-whitening filter that was applied
                %             To get mode2, filter the image thusly
                %               cMerge=ctf./pwf;
                %               filter=ctf1./cMerge;
                %               img=real(ifftn(fftn(orig_img).*ifftshift(filter)));
                %               ctf=cft1;
                %           Compute the mask
                msk=~meGetMask(mi,n0);  % Get the image mask
                %                 expandedMsk=BinaryConvolve(msk,circMask);
                expandedMsk=msk;
                
                %
                nPickerEntries=nPicks
                if totalNParts==0  % the first time, allocate part of the total stack.
                    totalStack=single(zeros(boxSize,boxSize,nPicks));
                    totalStackU=single(zeros(boxSize,boxSize,nPicks));
                end;
                nv=numel(mi.vesicle.x);
                vCoords=[mi.vesicle.x(:) mi.vesicle.y(:)];
                vRadii=mi.vesicle.r(:);
                stackSub=single(zeros(boxSize,boxSize,nPicks));
                stackImg=single(zeros(boxSize,boxSize,nPicks));
                sumImg=zeros(boxSize,boxSize);
                sumSub=zeros(boxSize,boxSize);
                disp('Extracting');
                %         figure(2);
                %         SetGrayscale;
                j=0;
                jPicks=0;
                nMasked=0;
                for i=1:nPicks  % scan the particle coordinates
                    type=mi.particle.picks(i,3);
                    vInd=mi.particle.picks(i,4);
                    skipVesicle=any(vInd==indBadVesicles);
                    %                     if skipVesicle
                    %                         type
                    %                     end
                    if (type>=16 && type<=32) && ~skipVesicle % valid particle flags
                        jPicks=jPicks+1;  % counter of all picked particles in micrograph
                        pCoords=mi.particle.picks(i,1:2);
                        iCoords=round(pCoords/ds)+1;
                        %                 Store for box display
                        iCoordsList(jPicks,:)=iCoords;
                        
                        j=j+1;  % counter of good particles
                        %                         if vInd<1  % no vesicle marked.  Find the nearest vesicle
                        %                             dists=sqrt(sum((vCoords-repmat(pCoords,nv,1)).^2,2));
                        %                             [minDist, vInd]=min(dists);
                        %                             if minDist > vRadii(vInd) % outside of membrane, check distance from mbn.
                        %                                 distm=abs(dists-vRadii);
                        %                                 [minDist, vInd]=min(distm);
                        %                             end;
                        %                             mi.particle.picks(i,4)=vInd;  % insert our estimated vesicle index.
                        %                             disp(['Inserted vesicle index for particle, vesicle ' num2str([i vInd])]);
                        %                         end;
                        %         Up to this point, all coordinates are in original pixel size.
                        intCoords=round(pCoords/ds);  % downsampled coordinates
                        fraCoords=pCoords/ds-intCoords;
                        xm=ExtractImage(mMerge,intCoords+1,padSize);
                        xmv=ExtractImage(mvMerge,intCoords+1,padSize);
                        shifts=fraCoords;
                        ximg=real(ifftn(fftn(xm).*fmasko.*FourierShift(padSize,-shifts)));
                        xmves=real(ifftn(fftn(xmv).*fmasko.*FourierShift(padSize,-shifts)));
                        %%
                        %             merge the padded particle images, reverse the contrast.
                        ximg2=-ximg.*padMask;
                        xsub2=-xmves.*padMask;
                        %                        Compute the angle relative to the y axis
                        %                           positive alpha means particle lies ccw from vertical
                        if rsMode
                            pVec=pCoords-vCoords(vInd,:);
                            alpha=atan2(pVec(2),pVec(1))-pi/2;
                            %        alpha=0;
                            %               rotate them
                            pSub=Crop(grotate(xsub2,-alpha),boxSize);
                            pImg=Crop(grotate(ximg2,-alpha),boxSize);
                        else
                            alpha=0;
                            pSub=Crop(xsub2,boxSize);
                            pImg=Crop(ximg2,boxSize);
                        end;
                        stackSub(:,:,j)=pSub;
                        sumSub=sumSub+pSub;
                        stackImg(:,:,j)=pImg;
                        sumImg=sumImg+pImg;
                        
                        subplot(2,4,3);
                        imacs(pSub);
                        title(mi.baseFilename,'interpreter','none');
                        subplot(2,4,4);
                        imacs(pImg);
                        subplot(2,4,7)
                        imacs(sumSub);  % unrotated image
                        title(j);
                        subplot(2,4,8);
                        imacs(sumImg);
                        drawnow;
                        
                        jt=j+totalNParts;  % pointer to overall stack index
                        si.miIndex(jt)=     fileIndex;
                        si.miParticle(jt)=  i;
                        if rsMode
                            si.alpha0(jt)=      alpha*180/pi;
                            si.yClick(jt)=      hypot(pVec(1),pVec(2))/ds;  % in units of si.pixA
                            si.rVesicle(jt)=    mi.vesicle.r(vInd)/ds;
                            si.sVesicle(jt)=    mi.vesicle.s(vInd);
                        end;
                    end; % if type
                end; % for i
                nparts=j;
                npicks=jPicks;
                
                % %                 %         Display the picks
                % %                 %%
                % %                 drawColors=[0 1 .3; 1 1 .3];
                % %                 figure(2);
                % %                 subplot(1,1,1);
                % %                 cimg=zeros([n0 3],'uint8');
                % %                 grayImage=imscale(GaussFilt(mvMerge,.3),256,.001);
                % %                 maskColor=[0 .4 .3];  % subtractive red
                % %                 emColor=[.1 .2 .0];
                % %                 for iColor=1:3
                % %                     cimg(:,:,iColor)=rot90(uint8(grayImage.*(1-expandedMsk*emColor(iColor)).*(1-msk*maskColor(iColor))));
                % %                 end;
                % %                 image(cimg);
                % %                 hold on
                % %                 bs2=boxSize/2;
                % %                 for ip=1:nparts
                % %                     xs=iCoordsList(ip,1)+[-bs2 -bs2 bs2  bs2 -bs2];
                % %                     ys=n0(2)-iCoordsList(ip,2)+[-bs2  bs2 bs2 -bs2 -bs2];
                % %                     plot(xs,ys,'-','color',drawColors(badFlags(ip)+1,:));
                % %                 end;
                % %                 hold off
                % %                 axis off;
                % %                 nBad=sum(badFlags(1:nparts));
                % %                 titleString=sprintf('%s   %d good   %d bad\n',...
                % %                     mi.baseFilename,nparts,npicks-nparts);
                % %                 title(titleString,'interpreter','none');
                % %                 %%
                stackSub=stackSub(:,:,1:nparts);  % truncate the stack
                stackImg=stackImg(:,:,1:nparts);
                
                totalStack(:,:,totalNParts+1:totalNParts+nparts)=stackSub;
                totalStackU(:,:,totalNParts+1:totalNParts+nparts)=stackImg;
                totalNParts=totalNParts+nparts
                figure(1);
                subplot(2,4,7)
                imacs(sum(totalStack,3));  % unrotated image
                title(totalNParts);
                subplot(2,4,8);
                imacs(sum(totalStackU,3));
                drawnow;
                
                if ~restoreFromSiFile && updateMiFile
                    WriteMiFile(mi,[infoPath miFiles{fileIndex}]);
                    disp([infoPath miFiles{fileIndex} ' updated']);
                end;
                
                
                if writeIndividualStacks
                    outname=[mi.stackPath mi.baseFilename 'st.mrc'];
                    WriteMRC(stackSub,pixA,outname);
                    outname=[mi.stackPath mi.baseFilename 'stu.mrc'];
                    WriteMRC(stackImg,pixA,outname);
                end;
            end; % if doExtractParticles
        else
            disp(' -no particle picks found');
        end;  % if np
    end;  % if isa struct
    fileIndex=fileIndex+1;
end; % for fileIndex
%
if totalNParts>0 || ~doExtractParticles
    
    si.pixA=pixA;
    si.fHighPass=fHP*pixA;
    mi1=si.mi{1};
    if isfield(mi1,'ppVars')
        si.mbnOffset=mi1.ppVars.membraneOffsetA/pixA;
    else
        si.mbnOffset=defaultMembraneOffsetA/pixA;
    end;
    si.weights=ones(size(mi.doses),'single');
    np=totalNParts;
    % truncate the si arrays.
    si.miIndex=     si.miIndex(1:np,1);
    si.miParticle=  si.miParticle(1:np,1);
    si.alpha0=      si.alpha0(1:np,1);
    si.yClick=      si.yClick(1:np,1);
    si.rVesicle=    si.rVesicle(1:np,1);
    si.sVesicle=    si.sVesicle(1:np,1);
    si.ctfs = ctfs;
    si.ctf1s=ctf1s;
    si.pwfs=pwfs;
    if restoreFromSiFile
        si.activeFlags=oldActiveFlags;
        si.activeFlagLog=oldActiveFlagLog;
    else
        si.activeFlags=true(np,1);
        si.activeFlagLog={[date '  StackExtractor2']};
    end;
    %
    %%
    if usePWFilter
        pwString='f';
    else
        pwString='';
    end;
    corrString='';
    if ~doExtractParticles
        corrString='noParticles'
    end;
    sizeString=sprintf('p%d',boxSize);
    p=strfind(mi.baseFilename,'_');
    if numel(p)>1
        baseName=mi.baseFilename(1:p(2)-1);
    else
        baseName=mi.baseFilename;
    end;
    outname=[mi.stackPath baseName weightString modeString pwString sizeString corrString];
    siName=[outname 'tsi.mat'];
    if strcmp(siName,oldSiName)
        ch=MyInput(['Overwrite the file ' siName],'n');
        if ch~='y'
            outname=MyInput('New base name',[outname 'z']);
            siName=[outname 'tsi.mat']
        end;
    end;
    % save the stackInfo structure
    save([outname 'tsi.mat'],'si');
    disp(['Wrote the total stack info ' outname '.tsi.mat']);
    
    if totalNParts>0 % Actually particles to store
        WriteMRC(totalStack,pixA,[outname stackSuffix]);
        disp(['Wrote the stack ' outname stackSuffix]);
        if writeUnsubtractedStack
            WriteMRC(totalStackU,pixA,[outname ustackSuffix]);
            disp(['Wrote the stack ' outname ustackSuffix]);
        end;
    else
        disp('No particle stacks written');
    end;
else
    disp('--nothing written.');
    disp(' ');
end;
