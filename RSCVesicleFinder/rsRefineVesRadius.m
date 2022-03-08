function [miNew, imgc, vesFit]=rsRefineVesRadius(img,oldVes,mask,mi,vIndex,effCTF,pars)
% Fit the radius (i.e. shape) to a single vesicle vIndex in the (downsampled, vesicle-subtracted, CTF-filtered)
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
% pars.displayOn


maxShiftA=200;    % max shift in A in fitting radius
maxShiftFrac=0.5; % max shift compared to radius
initStepA=30;  % initial step is ~ membrane width
convCriterion=1e-3; % std of simplex variations
minUnmaskedFraction=0.3; % Don't mask at all if less than this remains of a vesicle.
maskWeight=0.2;  % Scaling of old fit to replace the masked region.
termStepFraction=.8; % decrement of simplex steps per term
roundStepFraction=.8; % decrement of simplex steps per round
ndis=size(effCTF,1);
stepNTerms=2;


if pars.finalTerms<3 % We never consider the 2nd Fourier component
    pars.finalTerms=1;
end;

% Set up the number of iterations and free parameters for each round
totalPars=CountPars(pars.finalTerms);  % overall number of simplex parameters
% totalTerms=2 yields only 1 parameter as we ignore the 1st angular term.


roundNTerms=2:stepNTerms:pars.finalTerms+stepNTerms-1;
roundNTerms(end)=pars.finalTerms;
roundNTerms(1)=1; % we always skip term 2, i.e. jump from 1 to 3.

nRounds=numel(roundNTerms);  % number of simplex rounds
nIters=roundNTerms*pars.nIters; % number of simplex iters

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
maxShiftPix=(maxShiftA/mi.pixA+maxShiftFrac*vr(1))/ds; % max shift in our pixels

% Get the membrane cross-section density
pixA=mi.pixA*ds;
vd1=mi.vesicleModel;
vd=meDownsampleVesicleModel(vd1,ds)*pixA;  % for new fit

approxLoc=floor([vx vy]); % ones-based coordinate
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
vv=ves0(:); % our current model of the old vesicle
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

if pars.displayOn % Show the region to be fitted.
    clf;
    subplot(3,2,1);
    imags(imgc); % original image
    title(['Vesicle ' num2str(vIndex) ' / ' num2str(numel(mi.vesicle.x))]);
    drawnow;
    subplot(3,2,2); % Show the same, but with masking and presub
    imags(imgc0);
    title('pre-subtracted');
    drawnow;
%     subplot(3,2,3);  % in case we are displaying in LLSFit
end;

% Establish the initial Simplex steps
startingSteps=zeros(1,totalPars);
r1Step=initStepA/pixA;
startingSteps(1)=r1Step/5;
for i=2:2:totalPars
    startingSteps(i:i+1)=startingSteps(i-1)*termStepFraction;
end;
startingSteps(1)=r1Step;

% We update the shift values by locating the CC peak; the amplitude (and dc
% offset) by linear least squares; and the radius expansion by Simplex.

% Here is the fitting loop.  We do multiple rounds, increasing the number
% of radius terms with each round

p=struct;  % Struct of fitting parameters

p.pixA=pixA;
p.displayOn=pars.displayOn;
p.finalTerms=pars.finalTerms;
p.totalPars=totalPars;
p.imgc=imgc;
p.pixA=pixA;
p.t=[ctr ctr]+fracLoc; % translation
p.vd=vd;  % vesicle model
p.vr=vr;    % initial radius vector, updated as we fit
p.effCTF=effCTF;
p.mask1=mask1;
p.ccMask=ccMask;
p.logPScaled=logP*pars.hfVar/max(abs(vs),pars.hfVar);  % scale up to match ccf
p.rConstraints=pars.rConstraints;
p.convCriterion=convCriterion;
p.simpMask=false(1,p.totalPars);

%     ----------------rounds of simplex fitting ----------------
results=cell(pars.extraRound+1,5);
p.marker=char(0);

for iRound=1:nRounds+pars.extraRound % we'll do one round if nTerms(1)=0
    j=max(1,1+iRound-nRounds); % for remembering results
    jRound=min(nRounds,iRound);
    p.initPars=unpack(p.vr,pars.finalTerms);
    p.activeTerms=roundNTerms(jRound); % gradually increase this with rounds.
    nActivePars=CountPars(p.activeTerms);
    % Create the simplex mask and step values for setting the active parameters
    p.simpMask(1:nActivePars)=true;
    p.simpSteps=startingSteps; % true for round 1
    if jRound>1
        %         We scale down the steps of all previously-fit parameters
        p.simpSteps(2:2:oldActivePars)=startingSteps(2:2:oldActivePars)*roundStepFraction;
        p.simpSteps(3:2:oldActivePars)=startingSteps(3:2:oldActivePars)*roundStepFraction;
        p.simpSteps(oldActivePars+2:nActivePars)=r1Step/5; % Give a big boost to new pars
    end;
    if iRound>nRounds % bang at last step
%         as the Simplex starts with our previous best value, the error is
%         guaranteed to be lower in this 
        p.simpSteps=startingSteps(1)*[2 .5*ones(1,nActivePars-1)];
        p.marker='*';
    else
        oldActivePars=nActivePars;
    end;
    %         -----------Do the fit-----------
    [vr,s,t,err,vesFit]=SimplexFit(p,nIters(jRound));
    results(j,:)={vr s t err vesFit}; % store to find the best
    p.vr=vr;

    if p.displayOn % Show the final fit for this round
        subplot(3,2,3);
        imags(vesFit);
        title(num2str(p.simpMask));
        subplot(3,2,4);
        if jRound==nRounds  % flash the final fit
            imags(vesFit);
            pause(0.05);
        end;
        imags(p.imgc-vesFit);
        drawnow;
    end;

end; % for iRound

% Find the best extra round
nj=1+pars.extraRound;
errs=zeros(nj,1);
for j=1:nj
    errs(j)=results{j,4};
end;
[minErr,ij]=min(errs);

if p.displayOn && pars.extraRound
    marker=char(1,nj);
    marker(j)='*';
for j=1:nj
    subplot(6,4,16+j);
    imags(imgc-results{j,5});
    title([num2str(errs(j),3) marker(j)]);
end;
drawnow;
end;

vr=results{ij,1};
s=results{ij,2};
t=results{ij,3};
err=minErr;

% update the info structure
miNew=mi;
miNew.vesicle.r(vIndex,:)=0;
miNew.vesicle.r(vIndex,1:p.finalTerms)=vr*ds;   % complex radius
miNew.vesicle.s(vIndex,1)=s; % single amplitude value.

localXY=[approxLoc+t-ctr-1 1]'; % zero-based loc
globalXY=pars.M*localXY;
miNew.vesicle.x(vIndex)=globalXY(1); % zero-based position in orig image
miNew.vesicle.y(vIndex)=globalXY(2);
miNew.vesicle.s(vIndex,1,1)=s;
miNew.vesicle.err(vIndex,1)=err;

% abs(miNew.vesicle.r(12,:))

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

function sps=unpack(vr,nTerms)
% Get a vector of simplex pars from the vr terms.
if nargin<2
    nTerms=numel(vr);
end;
sps=vr(1);  % initialize the simplex parameters (unpack them)
vr(nTerms+1)=0; % pad with zeros if needed.
% for radial term fitting:
% pars=[vr(1) re(vr(3)) im(vr(3)) re(vr(4)) im(vr(4))...]
% pars will have 1+max(0,2*(nTerms-2)) elements.
for j=3:nTerms
    sps(2*(j-2))=real(vr(j));
    sps(2*(j-2)+1)=imag(vr(j));
end;
end

function np=CountPars(nt) % Convert number of terms to number of pars
np=max(1,1+2*(nt-2));
end

% ------------------------------------------------------

function [newVr,s,t,lsErr,vFit]=SimplexFit(p,nIters)
%   Simplex fit of radius parameters, followed by amplitude fit.
%   Inside the simplex loop are calls to LLSFit
%   to tune the overall amplitude and the translation t.
%   Simplex is called with {p.initPars,p.simpSteps,p.simpMask} and vR is
%   updated with Pack(P1,p.nTerms).  
%   The squared residual
%   is returned along with the final filtered
%   vesicle model vFit.
boundsErr=0;

%     Simplex initialization for the round
P1=Simplex('init',p.initPars,p.simpSteps,p.simpMask);

%     Simplex loop
for iRad=1:nIters  % radius fitting
    p.vr=Pack(P1,p.finalTerms);
    if p.activeTerms>2 % with a non-round vesicle, give errors for terms too large.
        boundsErr=1e3*sum( abs(p.vr(3:end))/p.vr(1) > p.rConstraints(3:end) );
    end;
    [vFit,p]=LLSFit(p); % Get the amplitude and translation fit
    res=p.imgc(:)-vFit(:);
    lsErr=res'*res;
    P1=Simplex(boundsErr+lsErr);  % update the vesicle radius

    if (mod(iRad,10)==0)
        subplot(323);
        imags(vFit);
        title([num2str(p.simpMask) '  r1= ' num2str(P1(1)*p.pixA,3) ' Å']);
        subplot(324);
        imags(p.imgc-vFit);
        title([num2str(iRad) '   ' num2str(p.activeTerms) '/' num2str(p.finalTerms) ' terms' p.marker]);
        drawnow;
%         disp(P1);
    end;
    if Simplex('converged',p.convCriterion)
        break; % quit the for loop
    end;
end;

% if p.s>0
% % %     % Find the centroid and do one more linear fit
%     P1=Simplex('centroid');
%      p.vr=Pack(P1,p.finalTerms);
%     [vFit,p]=LLSFit(p); % Get the translation one last time
% end;

    s=p.s;
    t=p.t;
    newVr=p.vr;
end

% ------------------------------------------------------------

function [vFit,p]=LLSFit(p)
%    Do the translation and basic amplitude fits of the membrane model.
%    Add the displacement of the cc peak to p.t to give the new
%    translation p.t. Also compute a=[const vesAmp]' to give p.s,
%    the amplitude of the general,
%    filtered vesicle with (non-updated) parameters {p.vr,p.vd}.
%    Return the model with constant and amplitude applied (but translation
%    not applied.)
ndis=size(p.imgc,1);
ctr=floor(ndis/2+1);
% First, determine the origin
tShift=1;
iter=0;
while tShift>.1 && iter<5
    v=-VesicleFromModelGeneral(ndis,p.vr,p.vd,p.t);
    F=ones(ndis^2,2,'single');
    fvFilt=fftn(v).*ifftshift(p.effCTF);
    vFilt=real(ifftn(fvFilt)).*p.mask1; % we really should fft this....
        % Compute the CCF
    ccMask=circshift(p.ccMask,round(ctr-p.t));
    cc=fftshift(real(ifftn(fftn(p.imgc).*conj(fftn(vFilt))))).*ccMask;
    [~, xi, yi]=max2di(cc+p.logPScaled);
    tUpdate=[xi yi]-ctr;
    p.t=p.t+tUpdate;
    tShift=sqrt(tUpdate*tUpdate');
        iter=iter+1;
end;

%     Compute the 2-component least-squares to get the amplitude s
F=ones(ndis^2,2);
F(:,2)=vFilt(:);
a=lscov(double(F),double(p.imgc(:)));
    % p.a=LinLeastSquares(F,p.imgc(:));
p.s=a(2);
vFit=F*a;  % fitted vesicle (amplitude and const)
vFit=reshape(vFit,ndis,ndis);

end

