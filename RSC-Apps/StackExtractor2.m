% StackExtractor2
% Version 2 uses vesicle shifts in subtraction and particle extraction.
% For each info file, read the micrographs and weight them according to the
% merge coefficients. For each particle, extract its region from each
% micrograph, shifted according to the mi.vesicle.shiftX and shiftY values.
% and sum the extracted images.
% Rotate the merged particle to the standard position, if desired
% (writeIndividualStacks=1) write *st.mrc and *stu.mrc files (stack and
% unsubtracted stack) into the Stack/ directory.  Also store *tstack.mrc
% and *tustack.mrc, the merged subtracted and unsubtracted stacks from all
% images, and *tsi.mat which contains the total si structure. Contrast is
% reversed so protein is white. The weights variable determines which
% exposures are used.  If the weights are not all 1s, e.g. [0 1 0] for the
% second exposure only, the saved files have the form *w010st.mrc and
% *w010tsi.mat.

boxSize=128;  % Size of boxes to be extracted from merged images.
% maskRadiusA=100;  % Size of disc surrounding the particle that is not
maskRadiusA=0;  % Size of disc surrounding the particle that may not
%  overlap with a masked region in the micrograph; otherwise the particle
%  is not extracted.
fHP=.001;    % highpass Gauss filter in A^-1
nZeros=1;    % number of zeros to include in high-defocus images
% touched by a mask.
ds=2;        % downsampling of boxed particles from original micrograph
weights=[1 1]  % which exposures are used.
mergeMode=1;  % do basic merging
usePWFilter=1
doExtractParticles=1;
% weights=[1 1 1]  % which exposures are used, and how they're weighted.
shiftMax=20; % maximum allowed vesicle shift (orig. pixels) between exposures
cpe0=16;     % default value, if not already in mi file
resetBasePath=1;

restoreFromSiFile=0;

writeIndividualStacks=0;
writeUnsubtractedStack=1;
useCorrection=0;  % extra vesicle residual subtraction
setOk4=1;  % force vesicle.ok(:,4) to 1.
stackSuffix='tstack.mrc';
ustackSuffix='tustack.mrc';
useShifts=0;  % Don't use shift data
defaultMembraneOffsetA=52;
skipOutliers=1;

stackDir='Stack'

% drawColors=[0 1 0; 1 .7 0];
dsScale=ds^2/4;  % fudge for vesicle model, don't really understand yet.
dds=2;       % further downsampling for micrograph display
vindex=0;    % force all vesicles to be modeled
padSize=NextNiceNumber(boxSize*1.5);
% Make the upsampled pad mask, and the particle fourier mask
padMask=fuzzymask(padSize,2,padSize*.48,padSize*.04);
fmasko=ifftshift(fuzzymask(padSize,2,.45*padSize,.05*padSize));

oldSiName='';
if restoreFromSiFile
    [oldSiName,pa]=uigetfile('*si.mat','Select si file to read');
    oldSi=load([AddSlash(pa) oldSiName]);
    miFiles=oldSi.si.mi;
    oldActiveFlags=oldSi.si.activeFlags;
    oldActiveFlagLog=oldSi.si.activeFlagLog;
    %     clear si;
else
    [miFiles, pa]=uigetfile('*mi.txt','Select mi files','multiselect','on');
    if ~iscell(miFiles)
        miFiles={miFiles};
    end;
end;

[rootPath, infoPath]=ParsePath(pa);
if isnumeric(rootPath)  % user clicked Cancel
    return
end;
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
        origWeights=mi.weights;
        mi.weights=weights;  % Modify the weights.
        
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
            disp(['Mask radius in A / pixels: ' num2str([maskRadiusA round(maskRadiusA/pixA)])]);
            disp(['   box size in A / pixels: ' num2str([round(boxSize*pixA) boxSize])]);
            circMask=fuzzymask(n0,2,maskRadiusA/pixA,2)>0;
        end;

        si.mi{fileIndex}=mi;  % store a copy of the micrograph info        
        if nPicks>0 % there are particles, insert filter functions

            if usePWFilter
                %     particle-sized noise-whitening filter
                Tic=meGetNoiseWhiteningFilter(mi,boxSize,ds,nZeros,fHP*pixA);  % for particle ctf
            else
                Tic=1;
            end;
            %             Set up CTF and other filters for this micrograph
            freqs=RadiusNorm(boxSize)/pixA;
            [pcoeffs, pctf, pdctfs]=meComputeMergeCoeffs2(freqs,mi,nZeros,mergeMode);
            
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
            
            if doExtractParticles  % go ahead and process micrographs
                if usePWFilter
                    [Ti,filterOk]=meGetNoiseWhiteningFilter(mi,n0,ds,nZeros,fHP*pixA);  % for micrograph
                    if ~filterOk
                        disp('No specimen-noise model, using simple highpass filter.');
                    end;
                else
                    Ti=1;
                end;
                disp('Compute image set');
                mSet=meMakeMergeImageSet(mi,mi.cpe,ds,Ti,nZeros,skipOutliers,mergeMode); % downsampling by ds
                mMerge=sum(mSet,3);  % merged image for checking
                nim=size(mSet,3);
                mi.weights=mi.weights(1:nim);
                imSet=find(mi.weights>0);  % The images we'll actually use.
                figure(1);
                clf;
                SetGrayscale;
                subplot(1,2,1);
                imacs(BinImage(mMerge,dds));
                title(mi.baseFilename,'interpreter','none');
                drawnow;
                
                disp('Compute vesicle models');
                if setOk4
                    mi.vesicle.ok(:,4)=true;  % Ignore the 4th column of this.
                end;
                doCTF=1;
                %         We fit only the vesicles that have particles
                %         vindex=find(hist(mi.particle.picks(:,4),1:nPicks)>0);
                fitM=size(mMerge);
                %         vSet=meMakeModelVesicleSet(mi,fitM,vindex,doCTF,Ti,nZeros)*dsScale;
                %           make vesicles with amplitude fitting
                vSet0=meMakeModelVesicleSet(mi,mMerge,vindex,doCTF,Ti,nZeros,mergeMode);
                vsCorr=zeros(size(vSet0),'single');
                mvSet0=mSet-vSet0;  % subtracted images
                if useCorrection
                    imacs(GaussFilt(sum(mvSet0,3),.2));
                    title('before correction');
                    drawnow;
                    %         Get the residuals, which will be subtracted from individual
                    %         exposures.  the residual computation is a linear operation.
                    for i=1:nim
                        vsCorr(:,:,i)=meGetVesicleResiduals(mi,sum(mvSet0,3));
                    end;
                    vSet=vSet0+vsCorr;
                else
                    vSet=vSet0;
                end;
                
                mvMerge=sum(mSet-vSet,3);  % final merged image
                %         imacs(BinImage(mvMerge,dds));
                imacs(GaussFilt(mvMerge,.2));
                title(mi.baseFilename,'interpreter','none');
                drawnow;
                % mtf=TestMarkParticles(mi,mMerge);
                
                % compute the CTF, including the prewhitening filter effects,
                % for individual particles
                freqs=RadiusNorm(boxSize)/pixA;
                [pcoeffs, pctf, pdctfs]=meComputeMergeCoeffs2(freqs,mi,nZeros,mergeMode);
                
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
                expandedMsk=BinaryConvolve(msk,circMask);
                
                %
                nPickerEntries=nPicks
                if totalNParts==0  % the first time, allocate part of the total stack.
                    totalStack=single(zeros(boxSize,boxSize,nPicks));
                    totalStackU=single(zeros(boxSize,boxSize,nPicks));
                end;
                nv=numel(mi.vesicle.x);
                vCoords=[mi.vesicle.x(:) mi.vesicle.y(:)];
                vRadii=mi.vesicle.r(:,1);
                stackSub=single(zeros(boxSize,boxSize,nPicks));
                stackImg=single(zeros(boxSize,boxSize,nPicks));
                sumImg=zeros(boxSize,boxSize);
                sumSub=zeros(boxSize,boxSize);
                if useShifts && isfield(mi.vesicle,'shiftX') && numel(mi.vesicle.shiftX)>0
                    shiftX=mi.vesicle.shiftX;
                    shiftY=mi.vesicle.shiftY;
                    shiftOk=mi.vesicle.shiftOk;
                else
                    disp('No shift values have been set.');
                    shiftOk=false(nv,nim);
                    shiftX=zeros(nv,nim);
                    shiftY=zeros(nv,nim);
                end;
                bad=(abs(shiftX)>shiftMax) | (abs(shiftY)>shiftMax) | ~shiftOk;
                shiftX(bad)=0;
                shiftY(bad)=0;
                disp('Extracting');
                %         figure(2);
                %         SetGrayscale;
                j=0;
                jPicks=0;
                nMasked=0;
                for i=1:nPicks  % scan the particle coordinates
                    type=mi.particle.picks(i,3);
                    if (type>=16 && type<=32)  % valid particle flags
                        jPicks=jPicks+1;  % counter of all picked particles in micrograph
                        pCoords=mi.particle.picks(i,1:2);
                        iCoords=round(pCoords/ds)+1;
                        badFlag=expandedMsk(iCoords(1),iCoords(2));  % a bad point
                        %                 Store for box display
                        iCoordsList(jPicks,:)=iCoords;
                        badFlags(jPicks)=badFlag;
                        
                        if ~badFlag
                            j=j+1;  % counter of good particles
                            
                            vInd=mi.particle.picks(i,4);
                            if vInd<1  % no vesicle marked.  Find the nearest vesicle
                                dists=sqrt(sum((vCoords-repmat(pCoords,nv,1)).^2,2));
                                [minDist, vInd]=min(dists);
                                if minDist > vRadii(vInd) % outside of membrane, check distance from mbn.
                                    distm=abs(dists-vRadii);
                                    [minDist, vInd]=min(distm);
                                end;
                                mi.particle.picks(i,4)=vInd;  % insert our estimated vesicle index.
                                disp(['Inserted vesicle index for particle, vesicle ' num2str([i vInd])]);
                            end;
                            
                            xves=zeros(padSize,padSize,nim);
                            ximg=zeros(padSize,padSize,nim);
                            %         Up to this point, all coordinates are in original pixel size.
                            for k=imSet  % loop over active images
                                intCoords=round(pCoords/ds);  % downsampled coordinates
                                fraCoords=pCoords/ds-intCoords;
                                xm=ExtractImage(mSet(:,:,k),intCoords+1,padSize);
                                xv=ExtractImage(vSet(:,:,k),intCoords+1,padSize);
                                shifts=[shiftX(vInd,k) shiftY(vInd,k)]/ds+fraCoords;
                                ximg(:,:,k)=real(ifftn(fftn(xm).*fmasko.*FourierShift(padSize,-shifts)));
                                xves(:,:,k)=real(ifftn(fftn(xv).*fmasko.*FourierShift(padSize,-shifts)));
                            end;
                            %%
                            %             merge the padded particle images, reverse the contrast.
                            ximg2=-sum(ximg,3).*padMask;
                            xsub2=-sum(ximg-xves,3).*padMask;
                            %                        Compute the angle relative to the y axis
                            %                           positive alpha means particle lies ccw from vertical
                            pVec=pCoords-vCoords(vInd,:);
                            alpha=atan2(pVec(2),pVec(1))-pi/2;
                            %        alpha=0;
                            %               rotate them
                            pSub=Crop(grotate(xsub2,-alpha),boxSize);
                            stackSub(:,:,j)=pSub;
                            sumSub=sumSub+pSub;
                            
                            pImg=Crop(grotate(ximg2,-alpha),boxSize);
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
                            si.alpha0(jt)=      alpha*180/pi;
                            si.yClick(jt)=      hypot(pVec(1),pVec(2))/ds;  % in units of si.pixA
                            si.rVesicle(jt)=    mi.vesicle.r(vInd)/ds;
                            si.sVesicle(jt)=    mi.vesicle.s(vInd);
                        end; % if ~badFlag
                    end; % if type
                end; % for i
                nparts=j;
                npicks=jPicks;
                
                %         Display the picks
                %%
                drawColors=[0 1 .3; 1 1 .3];
                figure(2);
                subplot(1,1,1);
                cimg=zeros([n0 3],'uint8');
                grayImage=imscale(GaussFilt(mvMerge,.3),256,.001);
                maskColor=[0 .4 .3];  % subtractive red
                emColor=[.1 .2 .0];
                for iColor=1:3
                    cimg(:,:,iColor)=rot90(uint8(grayImage.*(1-expandedMsk*emColor(iColor)).*(1-msk*maskColor(iColor))));
                end;
                image(cimg);
                hold on
                bs2=boxSize/2;
                for ip=1:nparts
                    xs=iCoordsList(ip,1)+[-bs2 -bs2 bs2  bs2 -bs2];
                    ys=n0(2)-iCoordsList(ip,2)+[-bs2  bs2 bs2 -bs2 -bs2];
                    plot(xs,ys,'-','color',drawColors(badFlags(ip)+1,:));
                end;
                hold off
                axis off;
                nBad=sum(badFlags(1:nparts));
                titleString=sprintf('%s   %d good   %d bad\n',...
                    mi.baseFilename,nparts,npicks-nparts);
                title(titleString,'interpreter','none');
                %%
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
                
                mi.weights=origWeights;   % Replace the original weights.
                if ~restoreFromSiFile
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
    si.weights=weights;
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
    if all(weights==1)
        weightString='';
    else
        format=['w%02d' repmat('-%02d',1,numel(weights)-1)];
        weightString=[sprintf(format,weights*10)];
    end;
    if mergeMode==1
        modeString='';
    else
        modeString=sprintf('m%1d',mergeMode);
    end;
    if usePWFilter
        pwString='f';
    else
        pwString='';
    end;
    if useCorrection
        corrString='c2';
    else
        corrString='v2';
    end;
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