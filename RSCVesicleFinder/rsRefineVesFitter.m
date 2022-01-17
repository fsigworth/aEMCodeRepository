function [miNew, imgc, vesFit]=rsRefineVesRadius(img,oldVes,mask,mi,vIndex,effCTF,pars,displayOn)
% Fit the radius (i.e. shape) to the vesicle vIndex in the (downsampled, vesicle-subtracted, CTF-filtered)
% micrograph img. To do the fitting we'll add back the model of
% vesicle(vIndex) if it was pre-subtracted.
% effCTF has zero frequency at the origin.
% A small image of size ndis=size(effCTF), imgc, is extracted from img and this is the portion
% that is fitted.  The corresponding fit is vesFit. The default for displayOn = 1.
% The starting vesicle info is taken from entry vindex in the mi.vesicle
% arrays, and the refined information is updated into the mi copy, miNew.
% miNew.vesicle.s(vIndex,1) is set to the overall amplitude of this vesicle
% as determined by least squares.
% Also miNew.vesicle.err(vIndex) is set, to give an error estimate.
%
% pars.preservedTerms is the number of radius terms to preserve from
%   previous fitting, in initializing the fitting. Can be inf.
% pars.finalTerms is the number of terms total to be fitted.
%   nTerms(:)>1 include higher-order radius terms
% The other fields of pars are
% pars.hfVar
% pars.rConstraints
% pars.nIters
% pars.M is the affine matrix to go from zero-based coordinates in img to
% the global zero-based coordinates in the original micrograph.

maxShiftA=200;    % max shift in A in fitting radius
maxShiftFrac=0.3;
initStepA=30;  % initial step is ~ membrane width
convCriterion=1e-3; % std of simplex variations
minUnmaskedFraction=0.3; % Don't mask at all if less than this remains of a vesicle.
maskWeight=0.2;  % Scaling of old fit to replace the masked region.
stepFraction=0.5; % decrement of simplex steps per round
ndis=size(effCTF,1);
termStep=2; % number of terms to add in each round


if nargin<8
    displayOn=1;
end;

% Set up the number of iterations and free parameters for each round
totalPars=CountPars(pars.finalTerms);  % overall number of simplex parameters
% totalTerms=2 yields only 1 parameter as we ignore the 1st angular term.

roundNTerms=1:termStep:pars.finalTerms;

nRounds=numel(roundNTerms);  % number of simplex rounds
nIters=(1:nRounds)*pars.nIters; % number of simplex iters

sigmaT=4;  % SD of prior for shifts
% Set the prior for translations
logP=-Radius(ndis).^2./(2*sigmaT^2);  % unscaled log prior

ds=pars.M(1,1); % Downsampling factor

% Get the vesicle coordinates in the downsampled image
origXY=[mi.vesicle.x(vIndex); mi.vesicle.y(vIndex); 1];
localXY=pars.M\origXY+1; % i.e. inv(M)*origXY; one-based

vx=localXY(1);
vy=localXY(2);
numPreservedTerms=min(pars.finalTerms,pars.preservedTerms);

% set up the starting parameters vr (radius vector) and vs (ampl scalar)
vr=mi.vesicle.r(vIndex,:)/ds;
vr(pars.finalTerms+1)=0; % Procrustean match of vector size.
vr(pars.finalTerms+1:end)=[];

vr(numPreservedTerms+1:end)=0; % zero out non-preserved terms
vs=mi.vesicle.s(vIndex,1,1); % use only constant term for amplitude
maxShiftPix=(maxShiftA/mi.pixA+maxShiftFrac*vr0(1))/ds; % max shift in our pixels

% Get the membrane cross-section density
pixA=mi.pixA*ds;
vd1=mi.vesicleModel;
vd=meDownsampleVesicleModel(vd1,ds)*pixA;  % for new fit

approxLoc=floor([vx vy]);
fracLoc=[vx vy]-approxLoc;
% approxLoc is the approximate position of the vesicle the whole image,
% one-based.
% Get a centered, cropped portion of the image, initial vesicles, mask
imgc0=ExtractImage(img,approxLoc,ndis);
ves0=ExtractImage(oldVes,approxLoc,ndis);
mask0=ExtractImage(single(mask),approxLoc,ndis);

% Make a mask for the ccf to prevent excessive shifts
ccmaskR=maxShiftPix+2;  % radius: add 2 pixels for offset errors
ccMask=fuzzymask(ndis,2,ccmaskR,ccmaskR/4);

ctr=floor(ndis/2+1);

% -----------Deal with masking----------
vv=ves0(:);
unmaskedFraction=(vv'*(vv.*mask0(:)))/(vv'*vv);
if isnan(unmaskedFraction)
    unmaskedFraction=1;
end;
mi.vesicle.af(vIndex,1)=unmaskedFraction;  % put in the active fraction
% Do no masking if the majority would be masked. Or maybe we should just
% skip this vesicle??
if unmaskedFraction<minUnmaskedFraction
    mask0=1;
end;
% Here we use the old fit to partly fill in the masked region
imgc=(imgc0+ves0).*mask0+(1-mask0).*ves0*maskWeight;  % put the old vesicle into the masked region.
mask1=mask0+(1-mask0)*maskWeight;  % use this mask for fitting

% -----------construct the image to be fitted------------
% We assume that img has been filtered equivalently to effCTF.  That is, if
% prewhitening has been done to img, it should also be included in effCTF.

if displayOn % Show the region to be fitted.
    subplot(2,2,1);
    imags(imgc);
    title(['Vesicle ' num2str(vIndex)]);

    subplot(2,2,2); % Show the same, but with masking
    imags(imgc0);
    title(unmaskedFraction);
    drawnow;
    subplot(2,2,3);  % in case we are displaying in LLSFit
end;

% Establish the initial Simplex steps
startingSteps=zeros(totalPars,1);
startingSteps(1)=initStepA/pixA;
for i=2:2:totalPars
    startingSteps(i:i+1)=startingSteps(i-1)*termStepFraction;
end;

% We update the shift values by locating the CC peak; the amplitude (and dc
% offset) by linear least squares; and the radius expansion by Simplex.

% Here is the fitting loop.  We do multiple rounds, increasing the number
% of radius terms with each round

p=struct;  % Struct of fitting parameters

p.imgc=imgc;
p.ringPeaks=mi.vesicle.extraPeaks/ds;
p.ringSigma=mi.vesicle.extraSD/ds;
p.pixA=pixA;
p.t=[ctr ctr]+fracLoc; % translation
p.vd=vd;  % vesicle model
p.vr=vr;    % initial radius vector, updated as we fit
p.effCTF=effCTF;
p.mask1=mask1;
p.ccMask=ccMask;
p.logPScaled=logP*pars.hfVar/max(abs(vs),pars.hfVar);  % scale up to match ccf
p.convCriterion=convCriterion;
p.displayOn=displayOn;

%     ----------------rounds of simplex fitting ----------------
for iRound=1:nRounds % we'll do one round if nTerms(1)=0
    p.initPars=unpack(p.vr);
    p.activeTerms=roundNTerms(iRound);
    nActivePars=CountPars(p.nActiveTerms);
    p.rConstraints=pars.rConstraints(1:p.activeTerms);
    p.rConstraints=p.rConstraints(:)';
    % Create the simplex mask and step values for setting the active parameters
    p.simpMask=false(1,totalPars);
    p.simpMask(1:nActivePars)=true;
    p.simpSteps=startingSteps; % true for round 1
    if iRound>1
%         We scale down the steps of all previously-fit parameters
        p.simpSteps(2:2:oldNParsRound)=startingSteps(2:2:oldnParsRound)*roundStepFraction;
        p.simpSteps(3:2:oldNParsRound)=startingSteps(3:2:oldnParsRound)*roundStepFraction;
    end;
    %         -----------Do the fit-----------
    [p.vr,s,t,err,vesFit]=SimplexFit(p,nIters(iRound));


%     rvals=p.vr;
%     save v2 p s vesFit rvals
%     disp('saved v2.mat'); %%%%%%%%%%%%%%%%%%
%     %         p.vr=vr1;  % update the starting radius
% 
    if p.displayOn
        subplot(2,2,3);
        imags(vesFit);
        title(num2str(p.simpMask));
        subplot(2,2,4);
        if iRound==nRounds  % flash the final fit
            imags(vesFit);
            pause(0.05);
        end;
        imags(p.imgc-vesFit);
        %     subplot(2,2,3);
        %     plot([sect(cc) sect(logPScaled)]);  % plot the CC and the prior
        drawnow;
    end;
end;


% update the info structure
miNew=mi;
sh=approxLoc-ctr; % zero-based shift to move to full image
miNew.vesicle.r(vIndex,:)=0;
miNew.vesicle.r(vIndex,1:finalTerms)=p.vr*ds;   % complex radius
localXY=[t+sh 1]'; % zero-based
globalXY=pars.M*localXY;
miNew.vesicle.x(vIndex)=globalXY(1); % zero-based position in image
miNew.vesicle.y(vIndex)=globalXY(2);
miNew.vesicle.s(vIndex,1,1)=s;
miNew.vesicle.err(vIndex,1)=err;

end % of main function




% ------------------------------------------------------------

function complexR=Pack(pars,nTerms)
% pack the simplex parameters into the vr format of radial terms.
if nTerms==1
    complexR=pars(1);
else
    complexR=zeros(1,nTerms(1));
    complexR(1)=pars(1);
    for it=3:nTerms(1)
        complexR(it)=pars(2*(it-2))+1i*pars(2*(it-2)+1);
    end;
end;
end

% ---------------------------------------------------

function pars=unpack(vr,nTerms)
pars=vr(1);  % initialize the simplex parameters (unpack them)
% for radial term fitting:
% pars=[vr(1) re(vr(3)) im(vr(3)) re(vr(4)) im(vr(4))...]
% pars will have 1+max(0,2*(nTerms-2)) elements.
for j=3:nTerms
    pars(2*(j-2))=real(vr(j));
    pars(2*(j-2)+1)=imag(vr(j));
end;
end

function np=CountPars(nt) % Convert number of terms to number of pars
    np=max(1,1+2*(nt-2));
end;
% ------------------------------------------------------

function [newVr,s,t,lsErr,vFit]=SimplexFit(p,nIters)
%   Simplex fit of radius parameters, followed by amplitude fit.
%   Inside the simplex loop are calls to LLSFit
%   to tune the overall amplitude and the translation t.
%   Simplex is called with {p.initPars,p.simpSteps,p.simpMask} and vR is
%   updated with Pack(P1,p.nTerms).  Finally, call
%   AmpFitNX to update the amplitude parameters s.  The squared residual
%   is returned alonsg with the final filtered
%   vesicle model vFit.
boundsErr=0;

%     Simplex initialization for the round
P1=Simplex('init',p.initPars,p.simpSteps,p.simpMask);
%     Simplex loop
for iRad=1:nIters  % radius fitting
    p.vr=Pack(P1,p.nTerms);
    %     if iRad==1
    %         disp(abs(p.vr));
    %     end;
    if p.nTerms(1)>2
        boundsErr=1e3*sum( abs(p.vr(3:end))/p.vr(1) > p.rConstraints(3:end) );
    end;
    [vfit,a,t]=LLSFit(p); % Get the amplitude and translation fit
    if any(isnan(vfit(:)))
        break % skip fitting if we got an NaN
    end;
    p.a=a;
    p.t=t;
    res=p.imgc(:)-vfit(:);
    lsErr=res'*res;
    %               Simplex update
    P1=Simplex(boundsErr+lsErr);  % update the vesicle radius
    if Simplex('converged',p.convCriterion)
        break;
    end;
end;

vesFit=vfit;
save v3 vesFit a t p
disp('v3.mat saved.') ;% ---these pars give the correct result!

if ~any(isnan(vfit(:)))
    % Find the centroid and do one more linear fit
    P1=Simplex('centroid');
    p.vr=Pack(P1,p.nTerms);
    [vFit,a,p.t,s]=LLSFit(p); % Get the translation one last time
    %     [vFit,s]=AmpFitNX(p);
    res=p.imgc(:)-vFit(:);
    lsErr=res'*res;
    t=p.t;
    newVr=p.vr;
else
    newVr=p.vr;
    [vFit,s]=AmpFitNX(p); % go ahead and get extra terms
    s=0*s;
    lsErr=p.imgc(:)'*p.imgc(:); % zero model fit.
    vFit=0*vFit;
end;
end

% ------------------------------------------------------------

function [vfit,a,t,s1]=LLSFit(p)
%    Do the translation and basic amplitude fits of the membrane model.
%    Add the displacement of the cc peak to p.t to give the new
%    translation t. Also compute a=[const vesAmp]'
%    the simple fit of a constant term and the amplitude of the general,
%    filtered vesicle with (non-updated) parameters {p.vr,p.vd,p.t}.
%    Return the model with constant and amplitude applied (but translation
%    not applied.)
ndis=size(p.imgc,1);
ctr=floor(ndis/2+1);

v=-VesicleFromModelGeneral(ndis,p.vr,p.vd,p.t);
F=ones(ndis^2,2,'single');
fvfilt=fftn(v).*ifftshift(p.effCTF);
vfilt=real(ifftn(fvfilt)).*p.mask0;
% Compute the CCF
cc=fftshift(real(ifftn(fftn(p.imgc).*conj(fvfilt)))).*p.ccMask;
[~, xi, yi]=max2di(cc+p.logPScaled);
t=p.t+[xi yi]-ctr;  % add the incremental translation

%     Compute the 2-component least-squares
F(:,2)=vfilt(:);
warning('off','MATLAB:singularMatrix');
a=LinLeastSquares(F,p.imgc(:));
warning('on','MATLAB:singularMatrix');
s1=a(1);
vfit=F*a;  % fitted vesicle (amplitude and const)
vfit=reshape(vfit,ndis,ndis);
if p.displayOn>1 % show each iteration
    subplot(2,2,3);
    imags(vfit);
    title(num2str(p.simpMask));
    subplot(2,2,4);
    imags(p.imgc-vfit);  % show the residual
    drawnow;
    %     imags(p.imgc-vfit);
    %     drawnow;
end;
end

% -----------------------------------------------------

% function [vfit,s,lsErr]=AmpFitNX(p)
% %         Linear fit of p.imgc with general model {p.vr,p.vd,p.t} as filtered
% %         with p.effCTF, to determine amplitudes s having theta dependence
% %           returns s with terms s(2...) complex-valued.  Also fits rings
% %           with additional LS fitted amplitudes.  The number of terms
% %           fitted nsTerms is set by p.nTerms(2) so that in the end s has the size
% %           nsTerms * (nx+1), where nx = numel(p.ringPeaks).
% nx=numel(p.ringPeaks);
% nt=nx+1;
% nsTerms=p.nTerms(2);
% ndis=size(p.imgc,1);
% vArray=zeros([ndis ndis nx+1],'single');
% vArray(:,:,1)=-VesicleFromModelGeneral(ndis,p.vr,p.vd,p.t);
% vAFilt=vArray;
% exPos=p.ringPeaks;
% exSigma=p.ringSigma;
% Ndis=ndis;
% if nt>1
%     vArray(:,:,2:nt)=-p.pixA*VesicleFromRings(ndis,exPos,exSigma,p.vr,p.t);
% end;
% F=zeros(ndis^2,(2*nsTerms-1)*nt,'single');
% 
% for k=1:nt
%     vAFilt(:,:,k)=real(ifftn(fftn(vArray(:,:,k)).*ifftshift(p.effCTF))).*p.mask0;
% end;
% 
% % Initialize the complex exponential
% [~,theta2D]=Radius(Ndis,p.t);
% theta=theta2D(:);
% w=exp(-1i*theta);  % complex exponential
% 
% % Get powers of the complex exponential
% %         vW=reshape(vAFilt,prod(imgSize),nx+1);
% %         fit will be of dc, 2*nsTerms-1 other terms.
% F(:,1)=ones(Ndis^2,1);  % constant term
% for l=1:nt
%     lOffset=(l-1)*(2*nsTerms-1);
%     vW=reshape(vAFilt(:,:,l),Ndis^2,1);
%     F(:,2+lOffset)=vW;  % const amplitude term
%     for k=2:nsTerms
%         vW=vW.*w;
%         F(:,2*k-1+lOffset)=real(vW);
%         F(:,2*k+lOffset)=imag(vW);
%     end;
% end;
% warning('off','MATLAB:singularMatrix');
% rs=LinLeastSquares(F,p.imgc(:));
% warning('on','MATLAB:singularMatrix');
% vfit=F*rs;  % fitted vesicle (amplitude and const)
% res=p.imgc(:)-vfit;
% lsErr=res'*res;
% vfit=reshape(vfit,Ndis,Ndis);
% 
% rs(1)=[];  % remove the dc term
% if nsTerms==1
%     s=rs;
% else
%     s=zeros(nsTerms*nt,1,'single');
%     for l=1:nt
%         lOffset=(l-1)*nsTerms;
%         lOffset2=(l-1)*(2*nsTerms-1);
%         s(1+lOffset)=rs(1+lOffset2);
%         s(2+lOffset:nsTerms+lOffset)...
%             =(rs(2+lOffset2:2:2*nsTerms-2+lOffset2))...
%             +1i*(rs(3+lOffset2:2:2*nsTerms-1+lOffset2));
%     end;
% end;
% s=reshape(s,nsTerms,nx+1);
% end

