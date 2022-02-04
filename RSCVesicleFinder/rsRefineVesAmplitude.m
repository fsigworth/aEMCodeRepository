function [miNew, imgc, vesFit]=rsRefineVesAmplitude(img,oldVes,mask,mi,vIndex,effCTF,pars)
% Fit the amplitude to a single vesicle vIndex in the (downsampled, vesicle-subtracted, CTF-filtered)
% micrograph img. To do the fitting we'll add back the model of
% vesicle(vIndex) if it was pre-subtracted.
% effCTF has zero frequency at the origin.
% A small image of size ndis=size(effCTF), imgc, is extracted from img and this is the portion
% that is fitted.  The corresponding fit is vesFit. The default for displayOn = 1.
% The starting vesicle info is taken from entry vindex in the mi.vesicle
% arrays, and the refined information is updated into the mi copy, miNew.
% miNew.vesicle.s(vIndex,pars.finalTerms) is set to the amplitude of this vesicle
% as determined by least squares.
% % Also miNew.vesicle.err(vIndex) is set, to give an error estimate.
%
% pars.finalTerms is the number of terms total to be fitted.
% The other fields of pars are
% pars.M is the affine matrix that maps zero-based coordinates in img to
% the global zero-based coordinates in mi.vesicle.x and y.

minUnmaskedFraction=0.3; % Don't mask at all if less than this remains of a vesicle.
maskWeight=0.2;  % Scaling of old fit to replace the masked region.
ndis=size(effCTF,1);


if pars.finalTerms<3 % We never consider the 2nd Fourier component
    pars.finalTerms=1;
end;

ds=pars.M(1,1); % Downsampling factor

% Get the vesicle coordinates in the downsampled image
origXY=[mi.vesicle.x(vIndex); mi.vesicle.y(vIndex); 1];
localXY=pars.M\origXY+1; % i.e. inv(M)*origXY; one-based

vx=localXY(1);
vy=localXY(2);
vr=mi.vesicle.r(vIndex,:)/ds;

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

ctr=floor(ndis/2+1);

% -----------Deal with masking----------
if pars.doPreSubtraction
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
else % no pre-subtraction and therefore no ves0.
    mask1=1;
imgc=imgc0;
end;


% -----------construct the image to be fitted------------
% We assume that img has been filtered equivalently to effCTF.  That is, if
% prewhitening has been done to img, it should also be included in effCTF.

if pars.displayOn % Show the region to be fitted.
    subplot(3,2,1);
    imags(imgc); % original image with vesicle
    title(['Vesicle ' num2str(vIndex)]);
    drawnow;
    subplot(3,2,2); % fitting image, with presub
    imags(imgc0);
    title('pre-subtracted');
    drawnow;
    %     subplot(3,2,3);  % in case we are displaying in LLSFit
end;

p=struct;  % Struct of fitting parameters

% p.pixA=pixA;
% p.displayOn=pars.displayOn;
% p.finalTerms=pars.finalTerms;
% p.imgc=imgc;
% p.t=[ctr ctr]+fracLoc; % translation
% p.vd=vd;  % vesicle model
% p.vr=vr;    % radius vector
% p.effCTF=effCTF;
% p.mask1=mask1;

ringPeaks=pars.peakPositionA/pixA;
ringSigma=pars.peakSigmaA/pixA;

%     ----------------LS amplitude fitting ----------------
% [s,err,vesFit]=AmplitudeFit(p);

% function [s,lsErr,vFit]=AmplitudeFit(p)
nx=numel(pars.peakPositionA);
nt=nx+1;
nsTerms=pars.finalTerms;
vArray=zeros([ndis ndis nt],'single');
t=[ctr ctr]+fracLoc; % translation
vArray(:,:,1)=-VesicleFromModelGeneral(ndis,vr,vd,t);
vAFilt=vArray;
if nt>1 % introduce the ring vesicle fumctions
    vArray(:,:,2:nt)=-pixA*VesicleFromRings(ndis,ringPeaks,ringSigma,vr,t);
end;
nFitPars=(2*nsTerms-1)*nt;

for k=1:nt
    vAFilt(:,:,k)=real(ifftn(fftn(vArray(:,:,k)).*ifftshift(effCTF))).*mask1;
end;

% Initialize the complex exponential
[~,theta2D]=Radius(ndis,t);
theta=theta2D(:);
w=exp(-1i*theta);  % complex exponential

% Get powers of the complex exponential
%         vW=reshape(vAFilt,prod(imgSize),nx+1);
%         fit will be of dc, 2*nsTerms-1 other terms.
F=ones(ndis^2,nFitPars);  % allocate, and set the constant term
for l=1:nt
    lOffset=(l-1)*(2*nsTerms-1);
    vW=reshape(vAFilt(:,:,l),ndis^2,1);
    F(:,2+lOffset)=vW;  % const amplitude term
    for k=2:nsTerms
        vW=vW.*w;
        F(:,2*k-1+lOffset)=real(vW);
        F(:,2*k+lOffset)=imag(vW);
    end;
end;

%% ------linear least-squares------
rs=lscov(double(F),double(imgc(:)));
vFit=F*rs;  % fitted vesicle (amplitude and const)
res=imgc(:)-vFit;
lsErr=res'*res;
vesFit=reshape(vFit,ndis,ndis);

rs(1)=[];  % remove the dc term

% Pack the fitted parameters into the complex ampltudes
if nsTerms==1
    s=rs;
else
    s=zeros(nsTerms*nt,1,'single');
    for l=1:nt
        lOffset=(l-1)*nsTerms;
        lOffset2=(l-1)*(2*nsTerms-1);
        s(1+lOffset)=rs(1+lOffset2);
        s(2+lOffset:nsTerms+lOffset)...
            =(rs(2+lOffset2:2:2*nsTerms-2+lOffset2))...
            +1i*(rs(3+lOffset2:2:2*nsTerms-1+lOffset2));
    end;
end;


if pars.displayOn % Show the final fit for this round
    subplot(3,2,3);
    imags(vesFit);
    subplot(3,2,4);
    imags(imgc-vesFit);
    drawnow;
end;

% update the info structure
miNew=mi;
s=reshape(s,[1 nsTerms nt]);
miNew.vesicle.s(vIndex,pars.finalTerms,nt)=0; % expand the array if needed
miNew.vesicle.s(vIndex,1:pars.finalTerms,:)=s; % amplitude values.


