function [miNew, imgc, vesFit]=rsQuickFitVesicleN(img,mask,mi,mi0,vIndex,effCTF,hfVar,nTerms,displayOn)
% function [miNew,imgc,vesFit]=rsQuickFitVesicleN(img,mask,mi,mi0,vIndex,effCTF,hfVar,nTerms,displayOn)
% Given a subtracted image, restore one vesicle using the model from mi0
% and the vesicle parameters from mi0. With one vesicle restored, fit by
% least-squares a vesicle model with fractional shifts.
% effCTF has zero frequency at the origin.
% A small image of size ndis=size(effCTF) is extracted from img and this is the portion
% that is fitted.  The default for displayOn = 1.
% The starting vesicle info is taken from entry vindex in the mi.vesicle
% arrays, and the refined information is updated into the mi copy, miNew.
% nTerms=0 fit amplitude only
% nTerms=1 fit position and amplitude (default)
% nTerms=2 fit radius too.
% nTerms>2 include higher-order radius terms
% The optional 2nd element of nTerms specifies the number of terms to use in
% the expansion of the amplitude vesicle.s(vIndex).  Default it is the same
% as nTerms(1).
minUnmaskedFraction=0.3; % Don't mask at all if less than this remains of a vesicle.
maxShiftA=100;    % max shift in Å
maxShiftFrac=0.5; % additional shift as a fraction of radius
initStepA=200;  % initial step is ~ membrane width
nRoundIters=200; % no of iterations to do for each step
parDiffCriterion=0;
xPeakPositionA=[-30 -20 20 30];  % location of extra peak relative to mbn center
% xPeakPositionA=[];
xPeakSigmaA=5;
convCriterion=1e-4;
maskWeight=0.2;  % Scaling of old fit to replace the masked region.

displayOn=1
% mi.vesicle.s=mi.vesicle.s(:,1:6); %%%%%%%%
% nTerms=3; %%%%%%%%%%%
nTerms

ndis=size(effCTF,1);
if nargin<9
    displayOn=1;
end;
if nargin<8
    nTerms=1;  % don't fit radius
end;

mask=single(mask);


if numel(nTerms)<2
    nTerms(2)=ceil(nTerms);
end;

% no of iterations to do for each step
% Add ntStep terms to r with each round.  E.g. ntStep=2 then
% nt1=3,5,7...nTerms
ntStep=3;  % increment number of radius terms with each round
% ntStep=1
stepFraction=0.25;
% stepFraction=1

nPars=max(0,1+(nTerms(1)-2)*2);  % overall number of simplex parameters
nRounds=max(0,ceil((nTerms(1)-1)/ntStep));  % number of simplex rounds

nIters=ones(1:nRounds)*nRoundIters;

sigmaT=2;  % SD of prior for shifts
n=size(img,1);
ds=mi.imageSize(1)/n;

% Get the vesicle coordinates relative to ms
vx=(mi.vesicle.x(vIndex))/ds+1;  % start with zero-based coordinates
vy=(mi.vesicle.y(vIndex))/ds+1;
vr0=mi.vesicle.r(vIndex,:)/ds;
vs=mi.vesicle.s(vIndex,1);

maxShiftPix=(maxShiftA/mi.pixA+maxShiftFrac*vr0(1))/ds; % max shift in our pixels

% Get the membrane cross-section density
pixA=mi.pixA*ds;
vd1=mi.vesicleModel;
vd=meDownsampleVesicleModel(vd1,ds)*pixA;  % for new fit

vdOrig=mi0.vesicleModel;
vd0=meDownsampleVesicleModel(vdOrig,ds)*mi0.pixA*ds;  % for old restoration

approxLoc=round([vx vy]);
% this is the approximate position of the vesicle relative to
% the center of the full-size image

% Get a centered, cropped portion of the image.
imgc0=ExtractImage(img,approxLoc,ndis);
mask0=ExtractImage(mask,approxLoc,ndis);
% Make a mask for the ccf to prevent excessive shifts
ccmaskR=maxShiftPix+2;  % radius: add 2 pixels for offset errors
ccmask=fuzzymask(ndis,2,ccmaskR,ccmaskR/4);

% -----------compute the original vesicle for adding back----------
v=-vs*VesicleFromModelGeneral(ndis,vr0,vd0,[vx vy]-approxLoc+ndis/2+1);
vv=v(:);
unmaskedFraction=(vv'*(vv.*mask0(:)))/(vv'*vv);
mi.vesicle.af(vIndex,1)=unmaskedFraction;  % active fraction


% Do no masking if the majority would be masked.
if unmaskedFraction<minUnmaskedFraction
    mask0=1;
end;
% here we try using the old fit to fill in the masked region
imgc0=imgc0.*mask0+(1-mask0).*v*maskWeight;  % put the old vesicle into the masked region.
mask1=mask0+(1-mask0)*maskWeight;

% -----------construct the image to be fitted------------
vFiltU=real(ifftn(fftn(v).*effCTF));
vFilt=vFiltU.*mask0;
imgc0=imgc0.*mask0+vFilt.*(1-mask0)*maskWeight;

% We assume that img has been filtered equivalently to effCTF.  That is, if
% prewhitening has been done to img, it should also be included in effCTF.
imgc=(imgc0+vFilt);  % add back the old vesicle

if displayOn
    subplot(2,2,1);
    imacs(imgc);
    title(vIndex(1));
    
    subplot(2,2,2);
    imacs(imgc0);
    title(unmaskedFraction);
    drawnow;
end;

logP=-Radius(ndis).^2./(2*sigmaT^2);  % unscaled log prior
ctr=floor(ndis/2+1);

nt1=nTerms(1); % number of terms to keep in the radius vector
if nt1<3
    nt1=1;
end;
if numel(vr0)<nt1
    vr0(nt1)=0;
    vr=vr0;
else
    vr=vr0(1:nt1);
end;
vr(2:end)=0;  % keep only the 1st term

% initial value for linear fit of amplitude
F=ones(ndis^2,2);  % we'll use F(:,1) for the constant term.
% F(:,2) for the vesicle itself.
P=vr(1);  % initialize the simplex parameters (unpack them)
% P=[r(1) re(r(3)) im(r(3)) re(r(4)) im(r(4))...]
for j=3:nTerms(1)
    P(2*(j-2))=real(vr(j));
    P(2*(j-2)+1)=imag(vr(j));
end;
%
sumT=[ctr ctr];  % accumulate the translations
parDiff=0;
    frac=.5;
    startingSteps(1)=initStepA/pixA;
    for i=2:2:numel(P)
        startingSteps(i:i+1)=startingSteps(1)*frac/i;
    end;
% We update the shift values by locating the CC peak; the amplitude (and dc
% offset) by linear least squares; and the radius expansion by Simplex.

subplot(2,2,3);  % in case we are displaying in LLSFit

P1=P;


% Here is the fitting loop.  We do multiple rounds, increasing the number
% of radius terms with each round
oldnParsRound=0;

simpMask=false;
res=imgc(:);

for iRound=1:nRounds
    % Create the simplex mask and step values for setting the active parameters
    nTermsRound=min(nTerms(1),1+ntStep*iRound);
    nParsRound=1+(nTermsRound-2)*2;  % number of active parameters in each round
    simpMask=false(1,nPars);
    simpMask(1:nParsRound)=true;
    simpSteps=startingSteps;
    simpSteps(1:oldnParsRound)=startingSteps(1:oldnParsRound)*stepFraction;
    oldnParsRound=nParsRound;
    
    %     Simplex initialization for the round
    P1=Simplex('init',P1,simpSteps,simpMask);
    %     Simplex loop
    oldParDiff=inf;
    for iRad=1:nIters(iRound)  % radius fitting
        vr1=Pack(P1);
        vfit=LLSFit(vr1); % Get the amplitude and translation fit
        res=imgc(:)-vfit;
        err=res'*res;
        %               Simplex update
        P1=Simplex(err);  % update the vesicle radius
        parDiff=min(oldParDiff,sum(abs(P1-oldP1))*pixA);
%         disp(parDiff)
        if parDiff<parDiffCriterion || Simplex('converged',convCriterion)
            break;
        end;
    end;
end;
% Find the centroid and do one more linear fit
P1=Simplex('centroid');
vr1=Pack(P1);
[vfit,a]=LLSFit(vr1); % Get the amplitude and translation fit

[vfit1,s1]=AmpFitN(vr,sumT,imgc,nTerms(2));

[vfit,s,extraS]=AmpFitNX(vr1,vd,sumT,imgc,nTerms(2));

% -------------Create extra components--------------
% nxc=numel(xPeakPositionA);  % number of extra components
% nvdx=2*ceil((max(xPeakPositionA)+4*xPeakSigmaA)/pixA)+1;
% cvdx=ceil(nvdx/2);
% vdx=zeros(nvdx,nxc);
% % df=.1;
% % means=nvd*[df:(1-2*df)/(nxc-1):1-df]';
% % sds=df/2;
% for i=1:numel(xPeakPositionA)
% vdx(:,i)=Gaussian(nvdx,1,xPeakSigmaA/pixA,cvdx+xPeakPositionA(i)/pixA);
% end;
% vdx=pixA*vdx;  % model 1V innter potential
% 
% [vfit1,s]=AmpFitN(vr1,sumT,imgc,nTerms(2));
% [vfit,s,extraS]=AmpFitNX(vr1,vd,vdx,sumT,imgc,nTerms(2));
vesFit=reshape(vfit,ndis,ndis);
diffIm=reshape(res,ndis,ndis);

% v=VesicleFromModelGeneral(ndis,s,vd,sumT,s);


if displayOn
    %     disp(xShift);
    subplot(2,2,4);
    imacs(diffIm);
    subplot(2,2,3);
    imacs(imgc-vfit1);
    %     subplot(2,2,3);
    %     plot([sect(cc) sect(logPScaled)]);  % plot the CC and the prior
    drawnow;
end;

% vrC=Pack(vr);

% update the info structure
miNew=mi;

% Expand the r and s fields if needed.  Otherwise clear out the unused
% terms.
if size(miNew.vesicle.r,2)<nt1
    miNew.vesicle.r(:,nt1)=0;
end;

if size(miNew.vesicle.s,2)<nTerms(2)
    miNew.vesicle.s(:,nTerms(2))=0;
end;

sh=approxLoc-ndis/2-1;

miNew.vesicle.r(vIndex,:)=0;
miNew.vesicle.r(vIndex,1:nt1)=vr1*ds;   % complex radius
%     miNew.vesicle.r(vindex,1)=vr*ds;  % radius isn't changed
miNew.vesicle.x(vIndex)=(sumT(1)+sh(1)-1)*ds;  % zero-based position in image
miNew.vesicle.y(vIndex)=(sumT(2)+sh(2)-1)*ds;
miNew.vesicle.s(vIndex,:)=0;
miNew.vesicle.s(vIndex,1:nTerms(2))=s(1:nTerms(2));


miNew.vesicle.extraPeaks=xPeakPositionA/mi.pixA;
miNew.vesicle.extraSD=xPeakSigmaA/mi.pixA;

nxc=numel(xPeakPositionA);
if ~isfield(miNew.vesicle,'extraS') || ndims(miNew.vesicle.extraS)<3
    miNew.vesicle.extraS=zeros(vIndex,nTerms(2),nxc, 'single');
elseif any(size(miNew.vesicle.extraS)<[vIndex nTerms(2) nxc]) % too small
    miNew.vesicle.s(vIndex,nTerms(2),nxc)=0;
end;
miNew.vesicle.extraS(vIndex,1:nTerms(2),1:nxc)=extraS;


    function complexR=Pack(pars)
        complexR=pars(1);
        for it=3:nTerms(1)
            complexR(it)=pars(2*(it-2))+1i*pars(2*(it-2)+1);
        end;
    end

    function [vfit,a]=LLSFit(vr)
        %         Do the translation and amplitude fits
        v=-VesicleFromModelGeneral(ndis,vr,vd,sumT);
        fvfilt=fftn(v).*effCTF;
        vfilt=real(ifftn(fvfilt)).*mask1;
        % Compute the CCF
        cc=fftshift(real(ifftn(fftn(imgc).*conj(fvfilt)))).*ccmask;
        logPScaled=logP*hfVar/vs;  % scale up to match ccf
        [~, xi, yi]=max2di(cc+logPScaled);
        sumT=sumT+[xi yi]-ctr;
        F(:,2)=vfilt(:);
        warning('off','MATLAB:singularMatrix');
        a=LinLeastSquares(F,imgc(:));
        warning('on','MATLAB:singularMatrix');
        vfit=F*a;  % fitted vesicle (amplitude and const)
if displayOn
    subplot(2,2,3); %%%%%%%%%%%%%%%%
            imags(reshape(vfit,ndis,ndis));
            title(num2str(simpMask));
            subplot(2,2,4);
            imags(imgc-reshape(vfit,ndis,ndis));
            title(parDiff);
            drawnow;
end;
    end

    function [vfit,s]=AmpFitN(vr,sumT,img,nsTerms)
%         Linear fit of amplitude s having theta dependence
%           returns s with terms s(2...) complex-valued.
        v=-VesicleFromModelGeneral(ndis,vr,vd,sumT);
        fvfilt=fftn(v).*effCTF;
        vfilt=real(ifftn(fvfilt)).*mask1;
        
        % Initialize the complex exponential
        [~,theta2D]=Radius(ndis,sumT);
        theta=theta2D(:);
        w=exp(-1i*theta);  % complex exponential
        w=w(:);
        wi=ones([ndis^2 nsTerms],'single');
        % Get powers of the complex exponential
        vW=vfilt(:);
        F=zeros(ndis^2,2*nsTerms);
        F(:,1)=1;
        F(:,2)=vW;
        for k=2:nsTerms
            vW=vW.*w;
            F(:,2*k-1)=real(vW);
            F(:,2*k)=imag(vW);
        end;
        warning('off','MATLAB:singularMatrix');
        s=LinLeastSquares(F,img(:));
        warning('on','MATLAB:singularMatrix');
        vfit=F*s;  % fitted vesicle (amplitude and const)
        subplot(2,2,3); %%%%%%%%%%%%%%%%
        imags(reshape(vfit,ndis,ndis));
        s(1)=[];  % remove the dc term
        for k=2:nsTerms % convert the later terms to complex
            s(k)=s(2*k-2)+1i*s(2*k-1);
        end;
        vfit=reshape(vfit,ndis,ndis);
    end
    


    function [vfit,s,extS]=AmpFitNX(vr,vd,sumT,img,nsTerms)
%         Linear fit of amplitude s having theta dependence
%           returns s with terms s(2...) complex-valued.  Also fits rings
%           with additional LS fitted amplitudes.
        nx=numel(xPeakPositionA);
        nt=nx+1;
        ndis=size(img,1);
        vArray=zeros([ndis ndis nx+1]);
        vAFilt=vArray;
        vArray(:,:,1)=-VesicleFromModelGeneral(ndis,vr,vd,sumT);
        exPos=xPeakPositionA/pixA;
        exSigma=xPeakSigmaA/pixA;
        vArray(:,:,2:nt)=-pixA*VesicleFromRings(ndis,exPos,exSigma,vr,sumT);
%         for k=1:nx
%             vArray(:,:,k+1)=-VesicleFromModelGeneral(imgSize,vr,vdx(:,k),sumT);
%         end;
        for k=1:nt
            vAFilt(:,:,k)=real(ifftn(fftn(vArray(:,:,k)).*effCTF)).*mask1;
        end;

        % Initialize the complex exponential
        [~,theta2D]=Radius(ndis,sumT);
        theta=theta2D(:);
        w=exp(-1i*theta);  % complex exponential
        
        % Get powers of the complex exponential
%         vW=reshape(vAFilt,prod(imgSize),nx+1);
%         fit will be of dc, 2*nsTerms-1 other terms.
        F=zeros(ndis^2,(2*nsTerms-1)*nt);
        F(:,1)=1;
        for l=1:nt
            lOffset=(l-1)*(2*nsTerms-1);
            vW=reshape(vAFilt(:,:,l),ndis^2,1);
            F(:,2+lOffset)=vW;
            for k=2:nsTerms
                vW=vW.*w;
                F(:,2*k-1+lOffset)=real(vW);
                F(:,2*k+lOffset)=imag(vW);
            end;
        end;
        warning('off','MATLAB:singularMatrix');
        rs=LinLeastSquares(F,img(:));
        warning('on','MATLAB:singularMatrix');
        vfit=F*rs;  % fitted vesicle (amplitude and const)
        vfit=reshape(vfit,ndis,ndis);
        subplot(2,2,3); %%%%%%%%%%%%%%%%
        imags(vfit);
        rs(1)=[];  % remove the dc term
        if nsTerms==1
            s=rs;
        else
        s=zeros(nsTerms*nt,1);
        for l=1:nt
            lOffset=(l-1)*nsTerms;
            lOffset2=(l-1)*(2*nsTerms-1);
            s(1+lOffset)=rs(1+lOffset2);
            s(2+lOffset:nsTerms+lOffset)...
                =(rs(2+lOffset2:2:2*nsTerms-2+lOffset2))...
                +1i*(rs(3+lOffset2:2:2*nsTerms-1+lOffset2));
        end;
        end;
        extS=reshape(s(nsTerms+1:end),nsTerms,nx);
        s=s(1:nsTerms);
    end
end
