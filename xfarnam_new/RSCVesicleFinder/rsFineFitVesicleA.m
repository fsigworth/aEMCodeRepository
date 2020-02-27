function [miNew diffIm vesFit]=rsQuickFitVesicle(img,mi,mi0,vindex,effCTF,hfVar,displayOn)
% function miNew=rsQuickFitVesicle(img,mi,mi0,vindex,effCTF,displayOn)
% Given a subtracted image with one vesicle restored, fit by
% least-squares a vesicle model with fractional shifts.  The
% starting vesicle info is taken from entry vindex in the mi.vesicle
% arrays, and the refined information is written back into this place.
% A small image of size ndis=size(effCTF) is extracted from img and this is the portion
% that is fitted.  The default for displayOn = 1.

ndis=size(effCTF,1);
if nargin<7
    displayOn=1;
end;

niters=2;
sigmaT=2;  % SD of prior for shifts
n=size(img,1);
ds=mi.imageSize(1)/n;

vx=(mi.vesicle.x(vindex))/ds+1;  % assume zero-based coordinates
vy=(mi.vesicle.y(vindex))/ds+1;
vr=mi.vesicle.r(vindex)/ds;
vs=mi.vesicle.s(vindex);

% Get the membrane cross-section density
vd1=mi.vesicleModel;
vd=meDownsampleVesicleModel(vd1,ds)*mi.pixA*ds;

vdOrig=mi0.vesicleModel;
vd0=meDownsampleVesicleModel(vdOrig,ds)*mi0.pixA*ds;

approxShift=round([vx vy]-n/2);
% this is the approximate position of the vesicle relative to
% the center of the full-size image

% Get a centered, cropped portion of the image.
imgc0=Crop(circshift(img,-approxShift),ndis);

% add back the original vesicle
v=-vs*VesicleFromModel(ndis,vr,vd0,[vx vy]-approxShift+(ndis-n)/2);
vfilt=real(ifftn(fftn(v).*effCTF));
% We assume that img has been filtered equivalently to effCTF.  That is, if
% prewhitening has been done to img, it should also be included in effCTF.
imgc=imgc0+vfilt;

if displayOn
    subplot(2,2,1);
    imacs(imgc);
    title(vindex(1));
    
    subplot(2,2,2);
    imacs(imgc0);
    drawnow;
end;


% initial value for location
const=ones(ndis,ndis)*100;
F(:,2)=const(:);
    % stepsizes=[.2 .2 .2];
% P=[vr ndis/2+1 ndis/2+1];
% 
% P=Simplex('init',P,stepsizes);

ctr=floor(ndis/2+1);
    v=-vs*VesicleFromModel(ndis,vr,vd,[ctr ctr]);
    fvfilt=fftn(v).*effCTF;
    vfilt=real(ifftn(fvfilt));
    eRef=vfilt(:)'*vfilt(:);  % approx. power in the reference.
    
logP=-Radius(ndis).^2./(2*sigmaT^2);  % unscaled log prior
% logP=0;
sumT=[ctr ctr];
xShift=zeros(1,niters);
for iter=1:niters
    % Construct the vesicle and interior functions
    cc=fftshift(real(ifftn(fftn(imgc).*conj(fvfilt))));
    logPScaled=logP*hfVar/vs;
    [mxc xi yi]=max2di(cc+logPScaled);
    sumT=sumT+[xi yi]-ctr;
    xShift(iter)=sumT(1);
    v=-VesicleFromModel(ndis,vr,vd,sumT);
    fvfilt=fftn(v).*effCTF;
    vfilt=real(ifftn(fvfilt));
    F(:,1)=vfilt(:);
    a=LinLeastSquares(F,imgc(:));
    vfit=F*a;

    
    diff=imgc(:)-vfit;
    err=diff'*diff;
    
    vesFit=reshape(vfit,ndis,ndis);
    diffIm=reshape(diff,ndis,ndis);
    %     subplot(2,2,3);
    %     imacs(vesfit);

end;
if displayOn
%     disp(xShift);
    subplot(2,2,4);
    imacs(diffIm);
    subplot(2,2,3);
%     imacs(vesFit);
    plot([sect(cc) sect(logPScaled)]);  % plot the CC and the prior
    drawnow;
end;

%% Construct the model image
% v=-VesicleFromModel(ndis,vr,vd,sumT);
% vfilt=real(ifftn(fftn(v).*effCTF));
% F=[vfilt(:) const(:)];
% a=LinLeastSquares(F,imgc(:));
% vfit=F*a;  % Ignore the constant term
% diff=imgc(:)-vfit;
% vesFit=reshape(vfit,ndis,ndis);
% diffIm=reshape(diff,ndis,ndis);  % residual image

% update the info structure
miNew=mi;
sh=approxShift-ndis/2+n/2;
% miNew.vesicle.r(vindex)=vr*ds;  % radius isn't changed
miNew.vesicle.x(vindex)=(sumT(1)+sh(1)-1)*ds;  % zero-based position in image
miNew.vesicle.y(vindex)=(sumT(2)+sh(2)-1)*ds;
miNew.vesicle.s(vindex)=a(1);



% disp([P(1) P(2) P(3) a']);  % show the main fit

end
