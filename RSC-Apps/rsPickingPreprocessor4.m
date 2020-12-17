function rsPickingPreprocessor4(fname,pars)
% rsPickingPreprocessor4
% This program does the correlations for SimpleRSPicker, writing a file
% *rscc.mat which contains the information.  It also places
% Some of the 2D arrays written out in the *rscc.mat file.
% mxCC is the computed amplitude of the putative particle at each pixel
% position.  The amplitude is relative to the
%  vesicle amplitude.  It is zero outside the vesicle region.
%  amxTemplInds gives the index of the reference giving that value. This is the
%  index k into eigenSet.vList(:,k) or, equivalently,
%  vList(:,hemi,gamma,angs) where angs are enumerated according to
%  rsListHemisphereAngles.  It is zero outside vesicle regions.
% mxVesInds gives the best-match or else the nearest vesicle to each pixel.
% mxRsos is 1 for a right-side-out particle at that position.  It
% is zero for inside-out and is also zero outside vesicle regions.
% mxVars is the local variance, computed as thesquare of the cc of
% the image with the first 1 or a few terms of the eigenimage expansion.

doBatchProcessing=(nargin>0);
if nargin<2
    pars=struct;
end;

defPars.overwrite=1; % overwrite existing *rscc.mat files
defPars.mapMode='Kv';
defPars.doNoiseWhitening=0; % if 0, use only a HP filter.
defPars.onlyChangedFlags=0; % do only mi files where the log shows flags were changed.

showTemplates=0;
defPars.readSubtractedImage=1;
pwFiltPars=[.002 0; .01 .5];  % generic pw filter parameters
defocusfHP=.001; % highpass in A^-1
phasefHP=.003;  % highpass to use if phase plate is in use.

% overrideDsm=2;  % downsample the merged image by this factor.
overrideDsm=0;
overrideB=40;

% These variables are copied to the mi file in rspLoadFiles.
localVarRadius=100;  % angstroms
maskPaddingA=40;     % extra space around outer radius, should be greater than maxBob in picker.
membraneThicknessA=60;

defPars.outputImageSize=960;  % size of output images.
defPars.outputImageSize=768;
defPars.mergedImageSuffix='s';
defPars.mapMode='Kv';
defPars.checkLog=0;

pars=SetOptionValues(defPars,pars);


nterms=26;
% nterms=51;
% nterms=1;
simulateImage=0;

doRelVarFigure=0;  % special figure of variance vs. number of terms
% nterms=40;       % use more terms for this figure.


ppVals=struct;
mapPath=AddSlash(fileparts(which('arGetRefVolumes'))); % Get our directory

disp(['Map: ' pars.mapMode]);

mapPath=AddSlash(fileparts(which('arGetRefVolumes'))); % Get our directory

switch pars.mapMode
    case 'Kv'
        ppVals.membraneOffsetA = 52;  % membrane is 52 � above particle center for RSO particle
        %    Hence a particle center > r-52� implies an inside-out particle.
        %     In turn, the center of beta subunits is about 50 � below the particle
        %     center
        %
        %     mapName='/Users/fred/Structures/kv1.2/KvMap.mrc';
        %     mapName='~/Structures/kv1.2/KvMap.mrc'; % weak membrane subtraction
        mapName=[mapPath '/KvMapMbnSub.mrc'];  % stronger membrane subtraction.
        dimShift=2;
        symmetry=4;
    case 'AMPAR'  % AMPAR
        ppVals.membraneOffsetA = -70;  % Membrane center is this distance from particle center
        mapName=[mapPath '3KG2map58.mrc'];
        dimShift=0;
        symmetry=2;
    case 'slo2'
        mapName='Slo2AMbnSub.mrc';
        mapName=[mapPath mapName];
        ppVals.membraneOffsetA=26;
        dimShift=2;
        symmetry=4;
   case 'bk'
        mapName='bk_TM0d50_d2d14f.mrc';
        mapName=[mapPath mapName];
        ppVals.membraneOffsetA=26;
        dimShift=0;
        symmetry=4;
        mapPadding=1.35;
    case 'atpsynth'
        mapName=[mapPath 'ATPSynthSub4A.mrc'];
        dimShift=0;
        symmetry=1;
        mapPadding=1.25;
        ppVals.membraneOffsetA=-76;
        
    otherwise
        warning('Unrecognized mapMode');
end;
if ~exist(mapName,'file')
    error(['The 3D map, full name ' mapName ' doesn''t exist.'])'
end
% new tighter values
ppVals.nAlpha=36;
ppVals.nBeta=12;  % even is best; number of betas per hemisphere.
ppVals.nGamma=24/symmetry;



% mapName='/Volumes/TetraData/Structures/AMPAR/3KG2mapsub5.8A.mrc';

% Have the user select some mi files: boilerplate
if ~exist('fname','var') || numel(fname)<1 || ~exist('doBatchProcessing','var') || ~doBatchProcessing
    [fname, pa]=uigetfile('*mi.txt','Select mi files','multiselect','on');
    if isnumeric(pa) % File selection cancelled
        return
    end;
    [rootPath, infoPath]=ParsePath(pa);
    if ~iscell(fname)
        fname={fname};
    end;
    cd(rootPath);
    doBatchProcessing=0;
else
    doBatchProcessing=1;
    rootPath=pwd;
end;

oldPixA=0;
figure(1);

%% ---------------Loop over micrographs--------------
for fileIndex=1:numel(fname) % Operate on a single micrograph
    tic
    if doBatchProcessing
      miName=fname{fileIndex};
    else
      miName=[infoPath fname{fileIndex}];
    end;
    disp([num2str(fileIndex) ' Reading ' miName]);
    mi=ReadMiFile(miName);
    
    if pars.onlyChangedFlags && ~strncmp(mi.log{end},'SimpleRSP:ChangedFlags',22)
        disp('--no vesicle changes, skipped.');
        continue;
    end;
    mi.basePath=rootPath;
    logIndices=miDecodeLog(mi);
    [origImg, mergeFullPath, mergedImgOk]=meReadMergedImage(mi,0,pars.mergedImageSuffix);
    ourProcPath=AddSlash(fileparts(mergeFullPath));
    outFileName=[ourProcPath mi.baseFilename 'rscc.mat'];
    if (exist(outFileName,'file') && logIndices(8)>logIndices(5)) && ~pars.overwrite
        disp('  A recent RSCC file exists.  Skipping.');
        continue;
    end;
    if pars.checkLog && (logIndices(6)<logIndices(5) || logIndices(5)<1)
        % no pw filter or no vesicle refinement
        disp('  Not ready for picking.');
        continue;
    end;
    
    numVesicles=0;
    if isfield(mi,'vesicle')
        numVesicles=numel(mi.vesicle.x);
    end;
    if numVesicles>0
        %     Pick up the merged image
        if ~mergedImgOk
            disp('No merged image.  Skipping.');
            continue;
        end;
        %%
        %     Possibly downsample it
        n1=size(origImg);
        ds0=ceil(mi.imageSize(1)/pars.outputImageSize);  % intended downsampling from merged image
        ds1=mi.imageSize(1)/n1(1);
        dsm=ds0/ds1;  % downsampling factor of merged image
        if overrideDsm
            dsm=overrideDsm;
        end;
        n=n1/dsm;                      % Output size
        if dsm>1  % further downsampling
            m0=DownsampleGeneral(origImg,round(n));
        else
            m0=origImg;
        end;
        n=size(m0);
        ctr=n/2+1;  % n must be even
        ds=mi.imageSize(1)/n(1);  % downsampling relative to original image
        pixA=mi.pixA*ds;
        disp(['Working image size is ' num2str(n)]);
        
        goodVesicle=all(mi.vesicle.ok(:,2:3),2); % in-range and refined
        badVesicle=~mi.vesicle.ok(:,2) & mi.vesicle.ok(:,3);  % found but out of range.
        
        if pars.readSubtractedImage
            [origSubImg, mergeFullPath, subImgOk]=meReadMergedImage(mi,0,['v' pars.mergedImageSuffix]);
            if ~subImgOk
                disp('No subtracted image. Skipping.');
                continue;
            end;
            if dsm~=1  % further downsampling
                m1=Downsample(origSubImg,n);
            else
                m1=origSubImg;
            end;
            mVesGood=m0-m1;
            mVesBad=0*mVesGood;
            
        else % compute vesicles instead
            
            
            %% Pick vesicles to work on
            %   definition: vesicles.ok(i,:) = [aVesicle inRange refined -- ]
            disp('Making model vesicles');
            mVesGood=meMakeModelVesicles(mi,n,find(goodVesicle));
            mVesBad=meMakeModelVesicles(mi,n,find(badVesicle));
            ves=mVesGood+mVesBad;
            %         ves=meMakeModelVesicles(mi,n,find(mi.vesicle.ok(:,3))); % everything that could be refined.
            m1=m0-ves;
        end;
        
        % Noise-whiten the image ----
        useNoiseWhitening=pars.doNoiseWhitening && (numel(mi.noiseModelPars)>0);
        if isfield(mi.ctf(1),'phi') && mi.ctf(1).phi>0
            fHP=phasefHP;
            disp('HP filter set for phase-plate data');
        else
            fHP=defocusfHP;
        end;
        if  useNoiseWhitening
            disp('Noise whitening');
            pwfImg=meGetNoiseWhiteningFilter(mi,n,0,1,fHP*pixA);
            if numel(pwfImg)<2
                disp('Using a generic filter');
                freqs=RadiusNorm(n)/pixA+1e-6;  % in inverse �, prevent divide by zero
                pwfImg=GaussHPKernel(n,fHP*pixA);
                for i=1:size(pwFiltPars,1)  % product of (gauss + const) terms.
                    f0=pwFiltPars(i,1);
                    a=pwFiltPars(i,2);
                    h=exp(-(f0./freqs).^2);
                    pwfImg=pwfImg.*(a+(1-a)*h);
                end;
            end;
        else
            disp('No specimen-noise whitening');
            pwfImg=GaussHPKernel(size(m1),fHP*pixA);
        end;
        m2=real(ifftn(fftn(m1).*ifftshift(pwfImg)));

%         vCorr=meGetVesicleResiduals(mi,m2);
%         m2=m2-vCorr;
        %     m=m1-ves;  % vesicle subtraction
        %  --------Start making the figure--------
%         figure(1);
        colormap jet(256);
        mysubplot(221);
        imaga(imscale(GaussFilt(m2,.3),256,.001));  % prewhitened image
        axis off;
        title(['Prewhitened: ' fname{fileIndex}],'interpreter','none');
        drawnow;
        % ----
        
        %%  Make templates
        if pixA~=oldPixA  % We haven't already made templates of the correct size
            oldPixA=pixA;
            % Load the 3D map
            disp(['Loading the 3D map ' mapName]);
            [origMap, mpixA]=ReadEMFile(mapName);
            nt1=size(origMap,1)*mpixA/pixA;  % final effective map size
            nt=ceil(nt1/8)*8;
            [map, finalmag]=DownsampleGeneral(origMap,nt,mpixA/pixA);
            map=map*pixA;  % approx amplitude correction (V-A scaling)
            map=shiftdim(map,dimShift);
            % magnifications=[mpixA/pixA finalmag]
            
            %% Create the list of angles for the templates
            membraneOffset=ppVals.membraneOffsetA/pixA;
            membraneOffsetA=ppVals.membraneOffsetA;
            gammaStep=360/(symmetry*ppVals.nGamma);
            
            %         hemiAngles run from alpha=[0..360) and beta=[0..90)
            [hemiAngles, angleInds]=rsListHemisphereAngles(ppVals.nAlpha, ppVals.nBeta);
            nHemiAngles=size(hemiAngles,1);
            nHemi=2;  % both hemispheres
            
            angleList=zeros(nHemi,ppVals.nGamma,nHemiAngles,3,'single');
            for j=1:ppVals.nGamma;
                gamma=(j-1)*gammaStep;
                for k=1:nHemiAngles
                    angleList(1,j,k,:)=[hemiAngles(k,:) gamma];
                    angleList(2,j,k,:)=[hemiAngles(k,1) 180-hemiAngles(k,2) gamma];
                    %                 angleList(2,j,k,:)=[[0 180]-hemiAngles(k,:) gamma];
                end;
            end;
            nAngles=numel(angleList)/3;
            
            % angle list is of size
            % (nHemi x nGamma x nHemiAngles, 3)  where nHemi=2 is the
            % number of hemispheres, and the hemisphere index is
            % fastest-varying.  Thus
            % the beta angles alternate such that, if betastep=1, they are
            % are (0, 180) for each gamma, then (1, 179) for each gamma, up to 89, 91.
            % that is, the same projected position is described twice.
            
            %% Make the templates
            
            disp(['Making ' num2str(nAngles) ' templates']);
            
            % allTemplates=rsMakeTemplatesQuick(angleList,map);
            allTemplates=rsMakeTemplates(reshape(angleList,nAngles,3),map);
            %             toc
            allTemplates=reshape(allTemplates,nt,nt,nHemi,ppVals.nGamma,nHemiAngles);
            
        end;
        
        %% Filter the templates according to the CTF and prewhitening filter
        [nt, nt, nHemi, ppVals.nGamma, nHemiAngles]=size(allTemplates);
        nAngles=nHemi*ppVals.nGamma*nHemiAngles;
        ne=NextNiceNumber(nt*1.3);  % increase the size to allow CTF rings
        mi1=mi;
        if overrideB
            for i=1:numel(mi.ctf)
                mi.ctf(i).B=overrideB;
            end;
        end;
        ctf=meGetEffectiveCTF(mi1,ne,ds);  % put in dqe, pw filter.
         if useNoiseWhitening % make a pre-whitening filter for the references
            pwfRef=meGetNoiseWhiteningFilter(mi1,ne,ds,1,fHP*pixA);

            if numel(pwfRef)<2
                disp('Using a generic filter');
                freqs=RadiusNorm(ne)/pixA+1e-6;  % in inverse �, prevent divide by zero
%                pwfRef=ones(ne,ne,'single');
                pwfRef=GaussHPKernel([ne ne],fHP*pixA);
                for i=1:size(pwFiltPars,1)  % product of (gauss + const) terms.
                    f0=pwFiltPars(i,1);
                    a=pwFiltPars(i,2);
                    h=exp(-(f0./freqs).^2);
                    pwfRef=pwfRef.*(a+(1-a)*h);
                end;
            end;
        else
%            pwfRef=ones(ne,ne,'single');
                pwfRef=GaussHPKernel([ne ne],fHP*pixA);            
        end;
        ctf=ctf.*pwfRef;
        
        % Pad the templates to avoid ctf artifacts
        xTemplates=Crop(reshape(allTemplates,nt,nt,nAngles),ne,1);
        nim=size(xTemplates,3);
        %     operate with the CTF and mask
        pwfImg=ifftshift(ctf);
        msk=fuzzymask(ne,2,0.45*ne,.1*ne);
        for i=1:nim
            xTemplates(:,:,i)=real(ifftn(fftn(xTemplates(:,:,i)).*pwfImg)).*msk;
        end;
        xTemplates=reshape(xTemplates,ne,ne,nHemi,ppVals.nGamma,nHemiAngles);
        ntstr=num2str(nt);
        nestr=num2str(ne);
        disp(['Templates expanded from ' ntstr 'x' ntstr ' to ' nestr 'x' nestr]);
        
        %% Make the eigenreferences
        eigenSet=rsMakeEigenrefs(xTemplates,nterms);
        % figure(1);
        % ImagicDisplay(eigenSet.imgs,2);
        %         %%
        %         %   Plot the band between min and max termVar as a gray region
        %         mxTermVar=max(eigenSet.termVar,[],2);
        %         mnTermVar=min(eigenSet.termVar,[],2);
        %         figure(2);
        %         SetGrayscale;
        %         xstep=.2;
        %         ystep=.01;
        %         [x,y]=ndgrid(0:xstep:nterms-xstep,0:ystep:1-ystep);
        %         pmask=255+0*x;
        %         for i=1:nterms
        %             pmask(round(x)==i & y>mnTermVar(i) & y<mxTermVar(i))=100;
        %         end;
        %         imac(x(:,1),y(1,:),pmask);
        %         set(gca,'fontsize',12);
        %         xlabel('Number of terms');
        %         ylabel('Relative power');
        %         %%
        
        if showTemplates
            %%  % Show the templates and the reconstructions
            %         in alternating order, in the stack rImgAlt.
            if showTemplatesAlternating
                timgs=reshape(eigenSet.imgs,ne*ne,nterms);
                nG=1;
                nH=2;
                nAngs=nHemiAngles;
                rImgAlt=single(zeros(ne,ne,2*nH,nG,nAngs));
                for k=1:nHemiAngles
                    for j=1:nG
                        for i=1:nH
                            %             rImg(:,:,2*i-1,j,k)=reshape(timgs*eigenSet.vList(:,i,j,k),ne,ne).*squeeze(eigenSet.ampList(i,j,k));
                            rImgAlt(:,:,2*i-1,j,k)=reshape(timgs*eigenSet.vList(:,i,j,k),ne,ne);
                            rImgAlt(:,:,2*i,j,k)=xTemplates(:,:,i,j,k);
                        end;
                    end;
                end;
                figure(5);
                ImagicDisplay2(rImgAlt,2);
                
                % % %         Create a set of references for Fig. 4 of the paper
                kImgs=single(zeros(ne,ne,nHemiAngles));
                qList=squeeze(eigenSet.vList(:,1,1,:));
                for k=1:nHemiAngles
                    kImgs(:,:,k)=rot90(reshape(timgs*qList(:,k),ne,ne),3);
                end;
                figure(1);
            else  % just show the reconstructed templates
                %%  % Make a complete set of reconstructions for comparison
                %   and make the average power spectrum
                nTotal=nHemiAngles*ppVals.nGamma*nHemi;
                timgs=reshape(eigenSet.imgs,ne*ne,nterms);
                rImg=single(zeros(ne,ne,nTotal));
                vL=reshape(eigenSet.vList,nterms,nTotal);
                sp=zeros(ne,ne);
                sp0=zeros(ne,ne);
                for i=1:nTotal
                    %     img=reshape(timgs*vL(:,i),ne,ne)*eigenSet.ampList(i);
                    img=reshape(timgs*vL(:,i),ne,ne);
                    sp0=sp0+abs(fftn(xTemplates(:,:,i))).^2;
                    sp=sp+abs(fftn(img)).^2;
                    rImg(:,:,i)=img;
                end;
                spr=Radial(fftshift(sp/nTotal));
                spr0=Radial(fftshift(sp0/nTotal));
                figure(2);
                mysubplot(2,1,1);
                semilogy([spr0 spr]);
                legend('Original','Expansion');
                title('Average reference power spectra');
                mysubplot(2,1,2);
                plot([cumsum(spr0) cumsum(spr)]);
                figure(5);  % show the reconstructed image
                ImagicDisplay2(rImg);
                figure(1);
            end;
        end;
        % mic=rsSortVesicles(mi);  % make a copy with the vesicles sorted by position.
        mic=mi;
        
%         mysubplot(222);
%         imaga(imscale(GaussFilt(m0,.3),256,.0005));
%         axis off;
%         drawnow;
%       -------display------
        mysubplot(221);
        imags(GaussFilt(origImg,1000/size(origImg,1)));
        title(miName,'interpreter','none');

        globalMask=meGetMask(mi,n);
        
        %%  ----------Evaluate the cc for each vesicle----------
        disp('Evaluating cross-correlations');
        disp([' using ' num2str(nterms) ' terms']);
        %         figure(3); clf; SetGrayscale;
        %         imac(imscale(GaussHP(GaussFilt(origImg,.1),.005),256,.003));
        
        %         %% Pick vesicles to work on
        %         %   definition: vesicles.ok(i,:) = [aVesicle inRange refined -- ]
        %         goodVesicle=all(mic.vesicle.ok(:,2:3),2); % in-range and refined
        %         if forceMedianAmplitude  % force amplitudes to the median
        %             medianS=median(mic.vesicle.s(:,1));
        %             nt=size(mic.vesicle.s,2);
        %             mic.vesicle.s(goodVesicle,:)=mic.vesicle.s(goodVesicle,:)...
        %                 .*medianS./repmat(mic.vesicle.s(goodVesicle,1),1,nt);
        %         end;
        %         %     allVesicles=mic.vesicle.ok(:,1);  % every vesicle that was found
        %         badVesicle=~mic.vesicle.ok(:,2) & mic.vesicle.ok(:,3);  % found but out of range.
        %         %     goodVesInds=find(okVesicles);
        %         %     badVesInds=find(allVesicles & ~okVesicles);  % exist, but not ok.
        
        gMaxVals=zeros(n,'single');
        gMaxValsU=zeros(n,'single');
        %     gMaxNCC=zeros(n,'single');
        mxTemplInds=uint16(zeros(n));
        mxVesInds=uint16(zeros(n));
        mxVars=zeros(n,'single');
        mxRsos=zeros(n,'single');
        mxDist=single(ones(n))*max(n);
        nulls=zeros(n,'single');
        badVesMask=false(n);
        msklRadius=round(localVarRadius/pixA);
        membraneThicknessPix=ceil(membraneThicknessA/pixA);
        nl=2*ceil(1.2*msklRadius);
        mskl=fuzzymask(nl,2,msklRadius,0.12*msklRadius);  % local mask
        npts=sum(mskl(:));
        
%         nves=numel(goodVesicle);
nves=numVesicles;
        disp(['Total vesicle entries: ' num2str(nves)]);
        maskPadding=ceil(maskPaddingA/pixA);
        
%             Correction for padded ms.mrc images
        if isfield(mi,'useMicrographCoords') && mi.useMicrographCoords
            cropOffset=floor((mi.padImageSize-mi.imageSize)/2);
        else
            cropOffset=[0 0];
        end;
        
        
        for i=1:numVesicles  % loop over all possible vesicles
            ctr=round(([mic.vesicle.x(i) mic.vesicle.y(i)]+cropOffset)/ds+1);  % shift to 1-based coords
            xDist=Radius(n,ctr)-mic.vesicle.r(i)/ds;  % Get the distance map, distance from mbn center
            qDist=xDist<mxDist;
            
            if goodVesicle(i) && mic.vesicle.r(i)>ds*membraneOffset  % We search in the vicinity of this one
                %                 maskRadii( outside-out, inside-out, overall boundary )
                os=sign(membraneOffset);
                maskRadii=mic.vesicle.r(i)/ds+membraneOffset*[-1 1 os]...
                    +maskPadding;
                %                     +maskPadding*[1 0 1];  % no padding for inside-out radius
                % %   ?ignore localVar calculation? maskRadii(3)=max(msklRadius,maskRadii(3)+max(-membraneOffset,0));
                
                % ---------- Get the single-vesicle correlation function here ------------
                [mxValsU, mxInds, mxNCC, mxVals, mxRso, localVar]=...
                    rsVesicleCorrelation6(-m2,mic,i,membraneOffset,...
                    maskRadii,angleInds,eigenSet);
                
                % Get the distance function for evaluating which is the closest
                nv=size(mxVals,1);
                
                % Do a local averaging of the squared CC in the whole vesicle.
                %         var=mxVals.^2;
                h=ifftshift(Crop(mskl,nv));  % average over the local mask
                filtVar=real(ifftn(fftn(localVar).*fftn(h))).*fuzzymask(nv,2,max(maskRadii(1:2)),.5);
                
                %     Pad all the quantities to a full-sized image
                xVals=ExtractImage(mxVals,ctr,n,1);
                xValsU=ExtractImage(mxValsU,ctr,n,1);
                xTemplInds=ExtractImage(mxInds,ctr,n,1);
                xVar=ExtractImage(filtVar,ctr,n,1);
                %         xVar=ExtractImage(localVar,ctr,n,1);
                %             xNCC=ExtractImage(mxNCC,ctr,n,1);
                xRso=ExtractImage(mxRso, ctr,n,1);
                
                %       Incorporate into the composite images
                qVar=xVar>mxVars;  % switch function for variances
                qAmp=xValsU>gMaxValsU; % switch function for particle amps
                
%            elseif badVesicle(i)
            else
                xVar=nulls;
                xVals=nulls;
                xValsU=nulls;
                xRso=nulls;
                xTemplInds=zeros(n,'uint16');
                qVar=qDist;  % switch function for variance is just the distance
                qAmp=false(n);
                badVesMask(xDist<membraneThicknessPix/2)=true;
            end;
            
            if goodVesicle(i) || badVesicle(i)  % Update the full-image statistics
                mxVars(qVar)=round(xVar(qVar));  % This is used for the variance criterion in SimpleRSPicker.
                mxDist(qDist)=round(xDist(qDist));  % global distance map
                gMaxVals(qAmp)=xVals(qAmp); % update the values depending on strongest signal
                gMaxValsU(qAmp)=xValsU(qAmp); % update the values depending on strongest signal
                mxVesInds(qAmp)=i;  % the related vesicle indices
                mxTemplInds(qAmp)=xTemplInds(qAmp);  % the template indices
                mxRsos(qAmp)=xRso(qAmp); % the right-side-out flags.
            end;
            
            if mod(i,10)==0
%                 mysubplot(223);
% % %                 imac(max((globalMask)*.015,gMaxVals.*globalMask)*jetScale);  % Show the CC function
%                 imacs(gMaxVals);  % Show the CC function
%                 axis off;
%                 title([num2str(i) ' / ' num2str(nves)]);
%                 drawnow;
%                 mysubplot(224);
%                 imags(mxTemplInds);
        mysubplot(222);
        imacs(gMaxValsU);
        axis off;
        title(['gMaxValsU ' num2str(nves) ' vesicles']);

        mysubplot(224);
%         imagesc(gMaxValsU);
%         title('gMaxValsU');
%         imags(mxTemplInds);
        imags(mxRsos);
        axis off;
        
        mysubplot(223);
        mxvi=mxVesInds;
        mxvi(1)=nves;  % force colormap
        imacs(mxvi);
        axis off;
        drawnow;
        
            end;
            % imacs(mxRsos);
            % drawnow;
        end;  % loop over vesicles
        % %%  Zero the cc of bad vesicles
        gMaxVals(badVesMask)=0;
        
% %         mysubplot(224);
% %         imacs(sqrt(mxVars));
% %         axis off;
%         mysubplot(222);
%         imagesc(gMaxVals);
%         title(['gMaxVals ' num2str(nves) ' vesicles']);
% 
%         mysubplot(224);
% %         imagesc(gMaxValsU);
% %         title('gMaxValsU');
%         imagesc(mxTemplInds);
%         
%         mysubplot(223);
%         imagesc(mxVesInds);
        
        %         imacs(min(1,gMaxVals));
%         imac(max((globalMask)*.015,gMaxVals.*globalMask)*jetScale);  % Show the CC function
        %     figure(4); SetGrayscale;
        %     imacs(gMaxNCC);
        %     title('gMaxNCC');
        
        mxCC=round(gMaxVals,3);  % round to 3 decimal places
        mxCCU=round(gMaxValsU,6);
        eigenImgs=eigenSet.imgs;
        vList=single(eigenSet.vList);
        %         Quantize the vesicle models
        mVesGood=round(mVesGood,3);
        mVesBad=round(mVesBad,3);
        
        partRadius=18; % ????
        logTxt=['rsPickingPreprocessor4 ' TimeStamp];
        save(outFileName,'mxCC','mxCCU','mxVars','mxVesInds',...
            'mxDist','mxTemplInds','mxRsos','partRadius', 'membraneOffsetA','ds',...
            'badVesMask','eigenImgs','vList','angleList','ppVals','fHP',...
            'pwfRef','pwfImg','mVesGood','mVesBad','logTxt','m0','m1');
        disp(['written: ' outFileName]);
    mi.log{end+1,:}=logTxt;
    WriteMiFile(mi,miName); % miName includes the infoPath

    else
        disp('No vesicles found.');
    end;
    
    toc
    disp(' ');
    % disp('pausing');
    % pause(15);
    % disp('continuing');
end;
% figure(1); clf;
