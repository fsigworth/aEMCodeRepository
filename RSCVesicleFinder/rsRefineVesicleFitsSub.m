function [miNew,vesFit]=rsRefineVesicleFitsSub(miOld,m,pars,displayOn)
% Do the vesicle refinement.  miOld is the original mi file; miNew is a
% copy with model parameters vesicle.extraPeaks and vesicle.extraSD set.

% pars.rTerms is a vector of radii in A, setting the number of terms in the fit.
%  e.g. rTerms=[90 120 160 200 inf] causes all vesicles with r<90A to be
%  fitted only with one term (spherical)

doTweakAmplitudes=0;
doPreSubtraction=pars.doPreSubtraction;

doFitFullSize=1;
doFitRadius=1;
switch pars.fitMode
    case 'LinOnly'
        doFitRadius=0;
    case 'RadiusOnly'
        doFitFullSize=0;
    case 'RadiusAndLin'
    otherwise
        warning(['Unrecognized fitMode: ' pars.fitMode]);
end;

defaultSValue=.0005*prod(miOld.imageSize)/numel(m);  % Value to assign if s(1) is otherwise 0 or negative.

maxVesiclesToFit=inf;
maxMaskLayers=2;   % Don't include any masking beyond merge and beam
useOkField=1;      % refine every vesicle for which ok is true.
doDownsampling=1;  % Downsample for speed
disA=1200;          % size of displayed/fitted window in angstroms
fHP=.003;          % Gauss high-pass filter for fitting
nZeros=1;          % number of zeros used in the merged CTF

%       Calculate variance from .30 to .45 x Nyquist; this is faster than using RadialPowerSpectrum:
n=size(m);
% pick the frequency range with an annulus in freq domain
annulus=fuzzymask(n,2,0.225*n,.05*n)-fuzzymask(n,2,0.15*n,.05*n);
spc=annulus.*fftshift(abs(fftn(m)).^2)/(n(1)*n(2));
hfVar0=sum(spc(:))/sum(annulus(:));

% % Handle cases where the ok field is old style
% if ~isfield(miOld.vesicle,'ok') || numel(miOld.vesicle.ok)<numel(miOld.vesicle.x) % no ok field at all
%     miOld.vesicle.ok=true(numel(miOld.vesicle.x),1);
% end;
nv=size(miOld.vesicle.x,1);

% Create the new mi structure
miNew=miOld;
% Set the extra peaks
miNew.vesicle.extraPeaks=pars.extraPeaks;
miNew.vesicle.extraSD=pars.extraSD;
% Insert the new 'active fraction' field.
if ~isfield(miNew.vesicle,'af') || numel(miNew.vesicle.af)~=nv  % active fraction (= unmasked)
    miNew.vesicle.af=ones(nv,1,'single');
end;

% Get image and pixel sizes ns, pixA
% ns will be the size of the downsampled image ms; pixA its pixel size.
n=size(m,1);
ds0=miOld.imageSize(1)/n;  % downsampling factor of m
pixA0=miOld.pixA*ds0;    % actual pixel size of m
if doDownsampling
    % further downsample the merged image to about 10A per pixel, yielding the image ms
    %     pars.targetPixA=10;  % maximum pixel size
    ns=NextNiceNumber(n*pixA0/pars.targetPixA,5,4);  % multiple of 4, largest factor 5.
    if ns<n
        disp(['Downsampling the micrograph to ' num2str(ns) ' pixels.']);
        ms=Downsample(m,ns);
    else
        ns=n;
        ms=m;
    end;
    ds=ds0*n/ns;  % downsampling factor of ms relative to original images.
    pixA=ds*miOld.pixA;  % pixA in the image ms.
else  % use the original merged image scale
    ds=ds0;
    pixA=pixA0;
    ns=n;
    ms=m;
end;
hfVar=hfVar0*ds0/ds;  % reduced hf spectral density after downsampling


%%  Get the original subtraction, and modify the amplitudes if necessary.
% miNew should the same as miOld except for the fields
% vesicle.extraPeaks and vesicle.extraSD.
%
% Find the vesicles to fit
nVesicles=sum(miOld.vesicle.ok(:,1));
vesList=find(miOld.vesicle.ok(:,1));
nVesicles=min(nVesicles,maxVesiclesToFit);
disp([num2str(nVesicles) ' vesicles to fit.']);

% Allocate the new miNew.vesicle r and s fields
maxR=max(miNew.vesicle.r(vesList,1));
maxNTerms=find(maxR<pars.rTerms/miNew.pixA,1);
nSTerms=maxNTerms;
np=numel(miNew.vesicle.extraPeaks);
nVes=size(miOld.vesicle.r,1);
% miNew.vesicle.r=zeros(nVes,maxNTerms,'single');
% miNew.vesicle.r(:,1)=miOld.vesicle.r(:,1);  % copy the basic radius
miNew.vesicle.r=miOld.vesicle.r;  % copy the basic radius **all terms**
miNew.vesicle.s=zeros(nVes,nSTerms,np+1,'single');
sVals=miOld.vesicle.s(:,:,1);  % Copy the basic amplitude
sVals(isnan(sVals))=0;
sVals(sVals<=0)=defaultSValue;
miNew.vesicle.s=sVals;


%   Compute the old subtraction (downsampled size)
nsPW=meGetNoiseWhiteningFilter(miOld,ns,ds,nZeros,fHP);
nsCTF=meGetEffectiveCTF(miOld,ns,ds);
msf=real(ifftn(fftn(ms).*ifftshift(nsPW)));  % High-pass filtered image
if doPreSubtraction
    vs=meMakeModelVesicles(miOld,ns,vesList,0,0); % no ctf or prewhitening
    vsf=real(ifftn(fftn(vs).*ifftshift(nsPW.*nsCTF)));  % hp filtered model
else
    vs=0;
    vsf=0;
end;

% Get the CTF information for the fitting regions
nds=NextNiceNumber(disA/pixA);  % size of display/fitting image
ndsPW=meGetNoiseWhiteningFilter(miOld,nds,ds,nZeros,fHP);
ndsCTF=meGetEffectiveCTF(miOld,nds,ds);

if doFitFullSize  % Get full-sized CTF and PW functions
    nPW=meGetNoiseWhiteningFilter(miOld,n,ds0,nZeros,fHP);
    nCTF=meGetEffectiveCTF(miOld,n,ds0);
    mf=real(ifftn(fftn(m).*ifftshift(nPW)));  % High-pass filtered image
    
    if displayOn
        figure(1);
        clf;
        imags(GaussFilt(mf,.1*ds0));
        title(['Prewhitened image ' miOld.baseFilename],'interpreter','none');
        drawnow;
    end;
    %     Make a full-sized subtraction
    if doPreSubtraction
        v=meMakeModelVesicles(miOld,n,vesList,0,0);
        vf=real(ifftn(fftn(v).*ifftshift(nPW.*nCTF)));
        if displayOn
            imags(GaussFilt(mf-vf,.1*ds0));
            title(['Preliminary subtraction ' miOld.baseFilename],'interpreter','none');
            drawnow;
        end;
    else
        vf=0;
    end;
    nd=NextNiceNumber(disA/pixA0);  % size of fit for full-size image
    ndPW=meGetNoiseWhiteningFilter(miOld,nd,ds0,nZeros,fHP);
    ndCTF=meGetEffectiveCTF(miOld,nd,ds0);
end;

%   If amplitude values are ridiculous, use linear least-squares to adjust the model scaling
if doTweakAmplitudes
    msAmpScale=(vsf(:)'*msf(:))/(vsf(:)'*vsf(:));
    if msAmpScale < 0.5 || msAmpScale > 1.5
        disp(['Amp scale for initial subtraction: ' num2str(msAmpScale)]);
        if msAmpScale>1e-3 % don't allow ridiculous values; it should be close to 1.
            vsf=vsf*msAmpScale;  % scale up the model
            vf=vf*msAmpScale;
            miOld.vesicle.s=miOld.vesicle.s*msAmpScale; % and the amplitude parameters
        else
            warning('Amp scale is too small, ignored.');
        end;
    end;
end;

%         Get the masks
layers=1:min(maxMaskLayers,numel(miOld.mask));
msmask=meGetMask(miOld,ns,layers);
mmask=meGetMask(miOld,n,layers);


%%  Actual fitting is done here
if displayOn
    figure(2);
end;
miNew.vesicle.ok(:,3)=true;  % we'll mark unfitted vesicles here.

if pars.listFits
    disp('  ind   1000s    r(Å)     ok   nTerms ------- 100s/s(1) -------------');
    %        1    2.779     205   1 1 1 1   4   24.80   23.71    0.00    0.00   0
end;
figure(2);

for j=1:nVesicles
    ind=vesList(j);
    ok=miOld.vesicle.ok(ind,1);  % The vesicle exists
    % ok=ok && any(ind==14);  %% just do index
    if useOkField && ~ok
        miNew.vesicle.s(ind,:)=0;
        continue
    end;
    
    %             Look up the number of radius terms to use.
    %             rTerms is something like [150 200 250 300 350 400 inf];
    %             r<150Å gets 1 term; <200 gets 2 terms; etc.
    nTerms=find(miNew.vesicle.r(ind,1)<pars.rTerms/miNew.pixA,1);
    % Set the number of amplitude terms.
    nTerms(2)=max(1,ceil(nTerms(1)*pars.fracAmpTerms));
    
    %             We'll do multiple rounds of fitting.  We first see how much
    %             we'll be expanding the number of terms from previously.
    origNTerms=min(pars.limitOrigNTerms,sum(miOld.vesicle.r(ind,:)~=0));
    nRounds=max(1,ceil((nTerms(1)-origNTerms+1)/3));
    
    %             Set up the number of terms we'll fit for each round
    nTR=zeros(nRounds,2); % no. terms in round for [r s]
    for i=1:nRounds
        nTR(i,1)=max(origNTerms,ceil(i*nTerms(1)/3)); % radius terms
        nTR(i,2)=min(nTerms(2),ceil(i*nTerms(2)/3));  % amp terms
    end;
    nTR(nRounds,:)=[nTerms(1) nTerms(2)]; % Fit all in the last round.
    %%%% rConstraints set here.
    rConstraints=ones(nTerms(1),1);
    rConstraints(2:nTerms)=0.4./((2:nTerms).^2)';
    
    %-------------------Basic fit------------------
    if doFitRadius % we're doing nonlinear fit
        if doPreSubtraction
            %                 First, compute the old model of the one vesicle in question
            vs1=meMakeModelVesicles(miOld,ns,ind,0,0); % no ctf or prewhitening
            vs1f=real(ifftn(fftn(vs1).*ifftshift(nsPW.*nsCTF)));  % filtered model
        else
            vsf=0;
            vs1f=0;
        end;
        %      --------------nonlinear fitting---------------
        [miNew,fitIms,vesFits]=rsQuickFitVesicle2(msf-vsf,vs1f,msmask,miNew,...
            ind,ndsCTF.*ndsPW,hfVar,nTR,rConstraints,displayOn);
        %             [miNew,fitIms,vesFits]=rsQuickFitVesicle2(msf-vsf,vs1f,msmask,miNew,...
        %                 ind,ndsCTF.*ndsPW,hfVar,nTR,displayOn); %%%%%%% 2nd round
        if displayOn
            subplot(2,2,4);
            imags(fitIms-vesFits);
            title(nTerms);
            drawnow;
        end;
    end;
    if doFitFullSize % do a linear fit of the vesicle.
        if doPreSubtraction
            v1=meMakeModelVesicles(miOld,n,ind,0,0);
            v1f=real(ifftn(fftn(v1).*ifftshift(nPW.*nCTF)));  % filtered model
        else
            v1=0;
            v1f=0;
        end;
        % if displayOn
        %     subplot(2,2,3);
        %             imags(v1);
        %             subplot(2,2,2);
        %             imags(fitIm-v1f);
        %             drawnow;
        % end;
        %                 Do a fit of only the amplitude terms
        [miNew,fitIm,vesFit]=rsQuickFitVesicle2(mf-vf,v1f,mmask,miNew,...
            ind,ndCTF.*ndPW,hfVar,[0 nTR(end,end)],rConstraints,displayOn & ~doFitRadius);  % no display
        
        if displayOn  % update the subtracted image
            subplot(2,2,4);
            imags(GaussFilt(fitIm-vesFit,pixA/20)); % 20 A filter
            drawnow;
        end;
    else
        vesFit=vesFits;
    end;
    
    if pars.listFits
        str=sprintf('%4d %8.3f  %6.1d  %2d%2d%2d%2d  %2d  %6.2f  %6.2f  %6.2f  %6.2f  ',...
            ind, 1000*miNew.vesicle.s(ind,1),...
            round(miNew.vesicle.r(ind,1)*miNew.pixA),...
            miNew.vesicle.ok(ind,:), nTerms,...
            100*abs(miNew.vesicle.s(ind,3:end,1))/miNew.vesicle.s(ind,1,1));
        disp(str);
    end;
end;

%         Disallow negative amplitudes, and nan values
miNew.vesicle.s(isnan(miNew.vesicle.s))=0;
badS=miNew.vesicle.s(:,1,1)<=0;
miNew.vesicle.s(badS,:,:)=0;
miNew.vesicle.ok(:,3)=miNew.vesicle.ok(:,1) & ~badS(:);  % unrefinable vesicles are marked 0

numberRefined=sum(miNew.vesicle.ok(:,3))
numberGood=sum(all(miNew.vesicle.ok(:,1:3),2))  %    exists, in range, refined
miNew.vesicle.refined=1;


