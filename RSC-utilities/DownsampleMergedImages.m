function DownsampleMergedImages(inNames,ipars)
if nargin<1
    inNames=[];
end;
if nargin<2
    ipars=struct;
end;
% Downsample each image in the current directory.

% suffixes={'.mrc' 'z.tif'};
defPars.suffixes={'m.mrc'};
defPars.inputDir='Merged/';
defPars.outputSuffix='ms.mrc';
defPars.outputDir='Merged/';
defPars.overwrite=0;

defPars.writeZTiff=0;
defPars.writeMRC=1;
defPars.outputImageSize=960;
defPars.doReplace=0;

pars=SetOptionValues(defPars,ipars);

names=cell(0,1);
outNames=cell(0,1);
numNames=0;

if numel(inNames)<1
    d=dir(pars.inputDir);
    numD=numel(d);
    for i=1:numD
        nm=d(i).name;
        [ok,outName]=TestSuffixes(nm,pars.suffixes,pars.outputSuffix);
        fullInName=[pars.inputDir nm];
        fullOutName=[pars.outputDir outName];
        if ok && (pars.overwrite || ~exist(fullOutName,'file'))
            numNames=numNames+1;
            names{numNames,1}=fullInName;
            outNames{numNames,1}=fullOutName;
        end;
    end;
else
    numNames=numel(inNames);
    for i=1:numNames
        mi=ReadMiFile(inNames{i});
        names{i}=[pars.inputDir mi.baseFilename pars.suffixes{1}];
        fullOutName=[pars.outputDir mi.baseFilename pars.outputSuffix];
        outNames{i,1}=fullOutName;
    end;
end;


disp([num2str(numNames) ' files to convert in ' pwd]);
szo=pars.outputImageSize;
disp(['Downsampling to ' num2str(szo)]);
CheckAndMakeDir(pars.outputDir,1);

for i=1:numNames
    name=names{i};
    outName=outNames{i};
    [m,pixA,ok]=ReadEMFile(name);
    ds=size(m,1)/szo(1);
    if ds>=1
        md=Downsample(m,szo);
        if pars.writeZTiff
            ztPars.snrRatio=200;  % target noise-to-noise ratio
            ztPars.lfCutoff=.2; % radius of frequency domain not fitted
            WriteZTiff(md,ds*pixA,outName,ztPars);
        elseif pars.writeMRC
            WriteMRC(md,ds*pixA,outName);
        end;
        disp([name ' -> ' outName]);
        if pars.doReplace
            str=['mv ' outName ' ' name];
            system(str);
        end;
    else
        disp([name ' is too small to downsample: ' num2str(szo*ds)]);
    end;
    
end;
end


function [ok,outName]=TestSuffixes(name,sufs,outSuf)
% Check if any of the suffixes exist for the given name.  If so, return
% the corresponding name.
nmLength=numel(name);
sufInd=0;
for j=1:numel(sufs)
    suf=sufs{j};
    start=nmLength-numel(suf)+1;
    if start>0
        if strcmpi(name(start:end),suf)
            sufInd=j;
            break;
        end;
    end;
end;
ok=sufInd>0;
if ok
    nsuf=numel(sufs{sufInd});
    outName=[name(1:end-nsuf) outSuf];
else
    outName='';
end;
end

%
%
% for j=1:numel(suffixes)
%
%
% flags=false(numD,1);
% for i=1:numD
%     ptrs=
%
%
%
% for fileIndex=1:numel(fname) % Operate on a single micrograph
%     tic
%     miName=[infoPath fname{fileIndex}];
%     disp(['Reading ' miName]);
%     mi=ReadMiFile(miName);
%     mi.basePath=rootPath;
%
%     if isfield(mi,'vesicle')&& numel(mi.vesicle.x)>0
%         %     Pick up the merged image
%         [origImg, mergeFullPath,mergedImageOk]=meReadMergedImage(mi);
%         if ~mergedImageOk
%             disp(['No merged image found: ' mergeFullPath]);
%         else  % got a merged image.
%             %%
%             %     Possibly downsample it
%             n1=size(origImg);
%             ds0=ceil(mi.imageSize(1)/outputImageSize);  % intended downsampling from merged image
%             ds1=mi.imageSize(1)/n1(1);
%             dsm=ds0/ds1;  % downsampling factor of merged image
%             if overrideDsm
%                 dsm=overrideDsm;
%             end;
%             n=n1/dsm;                      % Output size
%             if dsm~=1  % further downsampling
%                 m0=DownsampleGeneral(origImg,round(n));
%             else
%                 m0=origImg;
%             end;
%             n=size(m0);
%             ctr=n/2+1;  % n must be even
%             ds=mi.imageSize(1)/n(1);  % downsampling relative to original image
%             pixA=mi.pixA*ds;
%
%             goodVesicle=all(mi.vesicle.ok(:,1:3),2); % exists, in-range and refined
%             badVesicle= all(mi.vesicle.ok(:,1:3)==[1 0 1],2); % exists, not in-range
%
%             if readSubtractedImage
%                 [origSubImg, mergeFullPath,subImgOk]=meReadMergedImage(mi,0,'v');
%                 if ~subImgOk
%                     disp(['Subtracted image not found: ' mergeFullPath]);
%                 else
%                     if dsm~=1  % further downsampling
%                         m1=Downsample(origSubImg,n);
%                     else
%                         m1=origSubImg;
%                     end;
%                     mVesGood=m0-m1;
%                     mVesBad=0*mVesGood;
%                 end;
%             else % compute vesicles instead
%                 %% Pick vesicles to work on
%                 %   definition: vesicles.ok(i,:) = [aVesicle inRange refined -- ]
%                 disp('Making model vesicles');
%                 mVesGood=meMakeModelVesicles(mi,n,find(goodVesicle));
%                 mVesBad=meMakeModelVesicles(mi,n,find(badVesicle));
%                 ves=mVesGood+mVesBad;
%                 %         ves=meMakeModelVesicles(mi,n,find(mi.vesicle.ok(:,3))); % everything that could be refined.
%                 m1=m0-ves;
%                 subImgOk=true;
%             end;
%             if subImgOk  % go ahead and process.
%                 %         useNoiseWhitening=(numel(mi.noiseModelPars)>0);
%                 if  useNoiseWhitening
%                     disp('Noise whitening');
%                     H=meGetNoiseWhiteningFilter(mi,n,0,1,fHP*pixA);
%                     if numel(H)<2
%                         disp('Using a generic filter');
%                         freqs=RadiusNorm(n)/pixA+1e-6;  % in inverse �, prevent divide by zero
%                         H=ones(n,'single');
%                         for i=1:size(pwFiltPars,1)  % product of (gauss + const) terms.
%                             f0=pwFiltPars(i,1);
%                             a=pwFiltPars(i,2);
%                             h=exp(-(f0./freqs).^2);
%                             H=H.*(a+(1-a)*h);
%                         end;
%                     end;
%                     m2=real(ifftn(fftn(m1).*ifftshift(H)));
%                 else
%                     disp('No specimen-noise whitening');
%                     m2=m1;
%                     H=1;
%                 end;
%                 if useVesicleResiduals
%                     vCorr=meGetVesicleResiduals(mi,m2);
%                     m2=m2-vCorr;
%                 end;
%                 %     m=m1-ves;  % vesicle subtraction
%                 %
%                 figure(1);
%                 colormap jet(256);
%                 mysubplot(221);
%                 imaga(imscale(GaussFilt(m2,.3),256,.001));
%                 axis off;
%                 title(['Prewhitened: ' fname{fileIndex}],'interpreter','none');
%                 drawnow;
%                 %%
%                 if pixA~=oldPixA  % We haven't already made templates of the correct size
%                     oldPixA=pixA;
%                     % Load the 3D map
%                     disp(['Loading the 3D map ' mapName]);
%                     [origMap, mpixA]=ReadEMFile(mapName);
%                     expandN=NextNiceNumber(size(origMap,1)*mapPadding);
%                     expandMap=Crop(origMap,expandN);
%                     nt1=size(expandMap,1)*mpixA/pixA;  % final effective map size
%                     nt=ceil(nt1/8)*8;
%                     [map, finalmag]=DownsampleGeneral(origMap,nt,mpixA/pixA);
%                     map=map*pixA;  % approx amplitude correction (V-A scaling)
%                     map=shiftdim(map,dimShift);
%                     % magnifications=[mpixA/pixA finalmag]
%
%                     %% Create the list of angles for the templates
%                     membraneOffset=ppVals.membraneOffsetA/pixA;
%                     membraneOffsetA=ppVals.membraneOffsetA;
%                     gammaStep=360/(symmetry*ppVals.nGamma);
%
%                     %         hemiAngles run from alpha=[0..360) and beta=[0..90)
%                     [hemiAngles, angleInds]=rsListHemisphereAngles(ppVals.nAlpha, ppVals.nBeta);
%                     nHemiAngles=size(hemiAngles,1);
%                     nHemi=2;  % both hemispheres
%
%                     angleList=zeros(nHemi,ppVals.nGamma,nHemiAngles,3,'single');
%                     for j=1:ppVals.nGamma;
%                         gamma=(j-1)*gammaStep;
%                         for k=1:nHemiAngles
%                             angleList(1,j,k,:)=[hemiAngles(k,:) gamma];
%                             angleList(2,j,k,:)=[hemiAngles(k,1) 180-hemiAngles(k,2) gamma];
%                             %                 angleList(2,j,k,:)=[[0 180]-hemiAngles(k,:) gamma];
%                         end;
%                     end;
%                     nAngles=numel(angleList)/3;
%
%                     % angle list is of size
%                     % (nHemi x nGamma x nHemiAngles, 3)  where nHemi=2 is the
%                     % number of hemispheres, and the hemisphere index is
%                     % fastest-varying.  Thus
%                     % the beta angles alternate such that, if betastep=1, they are
%                     % are (0, 180) for each gamma, then (1, 179) for each gamma, up to 89, 91.
%                     % that is, the same projected position is described twice.
%
%                     %% Make the templates
%
%                     disp(['Making ' num2str(nAngles) ' templates']);
%
%                     % allTemplates=rsMakeTemplatesQuick(angleList,map);
%                     allTemplates=rsMakeTemplates(reshape(angleList,nAngles,3),map);
%                     %             toc
%                     allTemplates=reshape(allTemplates,nt,nt,nHemi,ppVals.nGamma,nHemiAngles);
%
%                 end;
%
%                 %% Filter the templates according to the CTF
%                 [nt, nt, nHemi, ppVals.nGamma, nHemiAngles]=size(allTemplates);
%                 nAngles=nHemi*ppVals.nGamma*nHemiAngles;
%                 ne=NextNiceNumber(nt*1.3);  % increase the size to allow CTF rings
%                 mi1=mi;
%                 if overrideB
%                     for i=1:numel(mi.ctf)
%                         mi.ctf(i).B=overrideB;
%                     end;
%                 end;
%                 ctf=meGetEffectiveCTF(mi1,ne,ds);  % put in dqe, pw filter.
%                 if useNoiseWhitening % make a pre-whitening filter for the references
%                     pwfRef=meGetNoiseWhiteningFilter(mi1,ne,ds);
%                     if numel(pwfRef)<2
%                         disp('Using a generic filter');
%                         freqs=RadiusNorm(ne)/pixA+1e-6;  % in inverse �, prevent divide by zero
%                         pwfRef=ones(ne,ne,'single');
%                         for i=1:size(pwFiltPars,1)  % product of (gauss + const) terms.
%                             f0=pwFiltPars(i,1);
%                             a=pwFiltPars(i,2);
%                             h=exp(-(f0./freqs).^2);
%                             pwfRef=pwfRef.*(a+(1-a)*h);
%                         end;
%                     end;
%                 else
%                     pwfRef=ones(ne,ne,'single');
%                 end;
%                 ctf=ctf.*pwfRef;
%
%                 % Pad the templates to avoid ctf artifacts
%                 xTemplates=Crop(reshape(allTemplates,nt,nt,nAngles),ne,1);
%                 nim=size(xTemplates,3);
%                 %     operate with the CTF and mask
%                 H=ifftshift(ctf);
%                 msk=fuzzymask(ne,2,0.45*ne,.1*ne);
%                 for i=1:nim
%                     xTemplates(:,:,i)=real(ifftn(fftn(xTemplates(:,:,i)).*H)).*msk;
%                 end;
%                 xTemplates=reshape(xTemplates,ne,ne,nHemi,ppVals.nGamma,nHemiAngles);
%                 ntstr=num2str(nt);
%                 nestr=num2str(ne);
%                 disp(['Templates expanded from ' ntstr 'x' ntstr ' to ' nestr 'x' nestr]);
%
%                 %% Make the eigenreferences
%                 eigenSet=rsMakeEigenrefs(xTemplates,nterms);
%                 % figure(1);
%                 % ImagicDisplay(eigenSet.imgs,2);
%                 %         %%
%                 %         %   Plot the band between min and max termVar as a gray region
%                 %         mxTermVar=max(eigenSet.termVar,[],2);
%                 %         mnTermVar=min(eigenSet.termVar,[],2);
%                 %         figure(2);
%                 %         SetGrayscale;
%                 %         xstep=.2;
%                 %         ystep=.01;
%                 %         [x,y]=ndgrid(0:xstep:nterms-xstep,0:ystep:1-ystep);
%                 %         pmask=255+0*x;
%                 %         for i=1:nterms
%                 %             pmask(round(x)==i & y>mnTermVar(i) & y<mxTermVar(i))=100;
%                 %         end;
%                 %         imac(x(:,1),y(1,:),pmask);
%                 %         set(gca,'fontsize',12);
%                 %         xlabel('Number of terms');
%                 %         ylabel('Relative power');
%                 %         %%
%
%                 if showTemplates
%                     %%  % Show the templates and the reconstructions
%                     %         in alternating order, in the stack rImgAlt.
%                     if showTemplatesAlternating
%                         timgs=reshape(eigenSet.imgs,ne*ne,nterms);
%                         nG=1;
%                         nH=2;
%                         nAngs=nHemiAngles;
%                         rImgAlt=single(zeros(ne,ne,2*nH,nG,nAngs));
%                         for k=1:nHemiAngles
%                             for j=1:nG
%                                 for i=1:nH
%                                     %             rImg(:,:,2*i-1,j,k)=reshape(timgs*eigenSet.vList(:,i,j,k),ne,ne).*squeeze(eigenSet.ampList(i,j,k));
%                                     rImgAlt(:,:,2*i-1,j,k)=reshape(timgs*eigenSet.vList(:,i,j,k),ne,ne);
%                                     rImgAlt(:,:,2*i,j,k)=xTemplates(:,:,i,j,k);
%                                 end;
%                             end;
%                         end;
%                         figure(1);
%                         ImagicDisplay2(rImgAlt,2);
%
%                         % % %         Create a set of references for Fig. 4 of the paper
%                         kImgs=single(zeros(ne,ne,nHemiAngles));
%                         qList=squeeze(eigenSet.vList(:,1,1,:));
%                         for k=1:nHemiAngles
%                             kImgs(:,:,k)=rot90(reshape(timgs*qList(:,k),ne,ne),3);
%                         end;
%
%                     else  % just show the reconstructed templates
%                         %%  % Make a complete set of reconstructions for comparison
%                         %   and make the average power spectrum
%                         nTotal=nHemiAngles*ppVals.nGamma*nHemi;
%                         timgs=reshape(eigenSet.imgs,ne*ne,nterms);
%                         rImg=single(zeros(ne,ne,nTotal));
%                         vL=reshape(eigenSet.vList,nterms,nTotal);
%                         sp=zeros(ne,ne);
%                         sp0=zeros(ne,ne);
%                         for i=1:nTotal
%                             %     img=reshape(timgs*vL(:,i),ne,ne)*eigenSet.ampList(i);
%                             img=reshape(timgs*vL(:,i),ne,ne);
%                             sp0=sp0+abs(fftn(xTemplates(:,:,i))).^2;
%                             sp=sp+abs(fftn(img)).^2;
%                             rImg(:,:,i)=img;
%                         end;
%                         spr=Radial(fftshift(sp/nTotal));
%                         spr0=Radial(fftshift(sp0/nTotal));
%                         figure(2);
%                         mysubplot(2,1,1);
%                         semilogy([spr0 spr]);
%                         legend('Original','Expansion');
%                         title('Average reference power spectra');
%                         mysubplot(2,1,2);
%                         plot([cumsum(spr0) cumsum(spr)]);
%                         figure(1);  % show the reconstructed image
%                         ImagicDisplay2(rImg);
%                     end;
%                 end;
%                 % mic=rsSortVesicles(mi);  % make a copy with the vesicles sorted by position.
%                 mic=mi;
%
%                 %         mysubplot(222);
%                 %         imaga(imscale(GaussFilt(m0,.3),256,.0005));
%                 %         axis off;
%                 %         drawnow;
%                 mysubplot(221);
%                 imags(GaussFilt(origImg,.1));
%                 title(miName,'interpreter','none');
%
%                 globalMask=meGetMask(mi,n);
%                 %%  ----------Evaluate the cc for each vesicle----------
%                 disp('Evaluating cross-correlations');
%                 disp([' using ' num2str(nterms) ' terms']);
%                 %         figure(3); clf; SetGrayscale;
%                 %         imac(imscale(GaussHP(GaussFilt(origImg,.1),.005),256,.003));
%
%                 %         %% Pick vesicles to work on
%                 %         %   definition: vesicles.ok(i,:) = [aVesicle inRange refined -- ]
%                 %         goodVesicle=all(mic.vesicle.ok(:,2:3),2); % in-range and refined
%                 %         if forceMedianAmplitude  % force amplitudes to the median
%                 %             medianS=median(mic.vesicle.s(:,1));
%                 %             nt=size(mic.vesicle.s,2);
%                 %             mic.vesicle.s(goodVesicle,:)=mic.vesicle.s(goodVesicle,:)...
%                 %                 .*medianS./repmat(mic.vesicle.s(goodVesicle,1),1,nt);
%                 %         end;
%                 %         %     allVesicles=mic.vesicle.ok(:,1);  % every vesicle that was found
%                 %         badVesicle=~mic.vesicle.ok(:,2) & mic.vesicle.ok(:,3);  % found but out of range.
%                 %         %     goodVesInds=find(okVesicles);
%                 %         %     badVesInds=find(allVesicles & ~okVesicles);  % exist, but not ok.
%
%                 gMaxVals=zeros(n,'single');
%                 gMaxValsU=zeros(n,'single');
%                 %     gMaxNCC=zeros(n,'single');
%                 mxTemplInds=uint16(zeros(n));
%                 mxVesInds=uint16(zeros(n));
%                 mxVars=zeros(n,'single');
%                 mxRsos=zeros(n,'single');
%                 mxDist=single(ones(n))*max(n);
%                 nulls=zeros(n,'single');
%                 badVesMask=false(n);
%                 msklRadius=round(localVarRadius/pixA);
%                 membraneThicknessPix=ceil(membraneThicknessA/pixA);
%                 nl=2*ceil(1.2*msklRadius);
%                 mskl=fuzzymask(nl,2,msklRadius,0.12*msklRadius);  % local mask
%                 npts=sum(mskl(:));
%
%                 nves=numel(goodVesicle);
%                 disp(['Total vesicle entries: ' num2str(nves)]);
%                 maskPadding=ceil(maskPaddingA/pixA);
%
%                 for i=1:nves  % loop over all possible vesicles
%                     ctr=round([mic.vesicle.x(i) mic.vesicle.y(i)]/ds+1);  % shift to 1-based coords
%                     xDist=Radius(n,ctr)-mic.vesicle.r(i)/ds;  % Get the distance map, distance from mbn center
%                     qDist=xDist<mxDist;
%
%                     %                 maskRadii( outside-out, inside-out, overall boundary )
%                     os=sign(membraneOffset);
%                     maskRadii=mic.vesicle.r(i)/ds+membraneOffset*[-1 1 os]...
%                         +maskPadding;
%                     %                     +maskPadding*[1 0 1];  % no padding for inside-out radius
%                     % %   ?ignore localVar calculation? maskRadii(3)=max(msklRadius,maskRadii(3)+max(-membraneOffset,0));
%
%                     if goodVesicle(i) && all(maskRadii>0)  % We search in the vicinity of this one
%
%                         % ---------- Get the single-vesicle correlation function here ------------
%                         [mxValsU, mxInds, mxNCC, mxVals, mxRso, localVar]=...
%                             rsVesicleCorrelation6(-m2,mic,i,membraneOffset,...
%                             maskRadii,angleInds,eigenSet);
%                         %                 mxValsU is the absolute amplitude, should be about .001
%                         %                 mxVals is the absolut amplitude, divided by the vesicle
%                         %                 amplitude s0.
%
%                         % Get the distance function for evaluating which is the closest
%                         nv=size(mxVals,1);
%
%                         % Do a local averaging of the squared CC in the whole vesicle.
%                         %         var=mxVals.^2;
%                         h=ifftshift(Crop(mskl,nv));  % average over the local mask
%                         filtVar=real(ifftn(fftn(localVar).*fftn(h))).*fuzzymask(nv,2,max(maskRadii(1:2)),.5);
%
%                         %     Pad all the quantities to a full-sized image
%                         xVals=ExtractImage(mxVals,ctr,n,1);
%                         xValsU=ExtractImage(mxValsU,ctr,n,1);
%                         xTemplInds=ExtractImage(mxInds,ctr,n,1);
%                         xVar=ExtractImage(filtVar,ctr,n,1);
%                         %         xVar=ExtractImage(localVar,ctr,n,1);
%                         %             xNCC=ExtractImage(mxNCC,ctr,n,1);
%                         xRso=ExtractImage(mxRso, ctr,n,1);
%
%                         %       Incorporate into the composite images
%                         qVar=xVar>mxVars;  % switch function for variances
%                         qAmp=xValsU>gMaxValsU; % switch function for particle amps
%
%                     elseif badVesicle(i)
%
%                         xVar=nulls;
%                         xVals=nulls;
%                         xValsV=nulls;
%                         xRso=nulls;
%                         xTemplInds=zeros(n,'uint16');
%                         qVar=qDist;  % switch function for variance is just the distance
%                         badVesMask(xDist<membraneThicknessPix/2)=true;
%                     end;
%
%                     if goodVesicle(i) || badVesicle(i)  % Update the full-image statistics
%                         mxVars(qVar)=round(xVar(qVar));  % This is used for the variance criterion in SimpleRSPicker.
%                         mxDist(qDist)=round(xDist(qDist));  % global distance map
%                         gMaxVals(qAmp)=xVals(qAmp); % update the values depending on strongest signal
%                         gMaxValsU(qAmp)=xValsU(qAmp); % update the values depending on strongest signal
%                         mxVesInds(qAmp)=i;  % the related vesicle indices
%                         mxTemplInds(qAmp)=xTemplInds(qAmp);  % the template indices
%                         mxRsos(qAmp)=xRso(qAmp); % the right-side-out flags.
%                     end;
%
%                     if mod(i,10)==0
%                         %                 mysubplot(223);
%                         % % %                 imac(max((globalMask)*.015,gMaxVals.*globalMask)*jetScale);  % Show the CC function
%                         %                 imacs(gMaxVals);  % Show the CC function
%                         %                 axis off;
%                         %                 title([num2str(i) ' / ' num2str(nves)]);
%                         %                 drawnow;
%                         %                 mysubplot(224);
%                         %                 imags(mxTemplInds);
%                         mysubplot(222);
%                         imacs(gMaxValsU);
%                         axis off;
%                         title(['gMaxValsU ' num2str(nves) ' vesicles']);
%
%                         mysubplot(224);
%                         %         imagesc(gMaxValsU);
%                         %         title('gMaxValsU');
%                         %         imags(mxTemplInds);
%                         imags(mxRsos);
%                         axis off;
%
%                         mysubplot(223);
%                         mxvi=mxVesInds;
%                         mxvi(1)=nves;  % force colormap
%                         imacs(mxvi);
%                         axis off;
%                         drawnow;
%
%                     end;
%                     % imacs(mxRsos);
%                     % drawnow;
%                 end;  % loop over vesicles
%                 % %%  Zero the cc of bad vesicles
%                 gMaxVals(badVesMask)=0;
%                 gMaxValsU(badVesMask)=0;
%                 % %         mysubplot(224);
%                 % %         imacs(sqrt(mxVars));
%                 % %         axis off;
%                 %         mysubplot(222);
%                 %         imagesc(gMaxVals);
%                 %         title(['gMaxVals ' num2str(nves) ' vesicles']);
%                 %
%                 %         mysubplot(224);
%                 % %         imagesc(gMaxValsU);
%                 % %         title('gMaxValsU');
%                 %         imagesc(mxTemplInds);
%                 %
%                 %         mysubplot(223);
%                 %         imagesc(mxVesInds);
%
%                 %         imacs(min(1,gMaxVals));
%                 %         imac(max((globalMask)*.015,gMaxVals.*globalMask)*jetScale);  % Show the CC function
%                 %     figure(4); SetGrayscale;
%                 %     imacs(gMaxNCC);
%                 %     title('gMaxNCC');
%
%                 mxCC=round(gMaxVals,3);  % round to 3 decimal places
%                 mxCCU=round(gMaxValsU,6);
%                 eigenImgs=eigenSet.imgs;
%                 vList=single(eigenSet.vList);
%                 %         Quantize the vesicle models
%                 mVesGood=round(mVesGood,3);
%                 mVesBad=round(mVesBad,3);
%                 log=['rsPickingPreprocessor4 ' TimeStamp];
%
%                 partRadius=18; % ????
%                 save([mi.procPath mi.baseFilename 'rscc.mat'],'mxCC','mxCCU','mxVars','mxVesInds',...
%                     'mxDist','mxTemplInds','mxRsos','partRadius', 'membraneOffsetA','ds',...
%                     'badVesMask','eigenImgs','vList','angleList','ppVals','pwfRef',...
%                     'mVesGood','mVesBad','log');
%                 disp(['written: ' mi.procPath mi.baseFilename 'rscc.mat']);
%             else
%                 disp('No vesicles found.');
%             end;
%             toc
%             disp(' ');
%             % disp('pausing');
%             % pause(15);
%             % disp('continuing');
%         end; % if subImgOk
%     end; % if mergeImageOk
% end;
% % figure(1); clf;
