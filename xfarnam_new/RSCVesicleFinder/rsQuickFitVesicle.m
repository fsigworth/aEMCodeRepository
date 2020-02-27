function [miNew, imgc, vesFit]=rsQuickFitVesicle(img,mask,mi,mi0,vindex,effCTF,hfVar,mode,displayOn)
% function [miNew,imgc,vesFit]=rsQuickFitVesicle(img,mask,mi,mi0,vindex,effCTF,hfVar,mode,displayOn)
% Given a subtracted image, restore one vesicle using the model from mi0
% and the vesicle parameters from mi. With one vesicle restored, fit by
% least-squares a vesicle model with fractional shifts.
% effCTF has zero frequency at the origin.
% A small image of size ndis=size(effCTF) is extracted from img and this is the portion
% that is fitted.  The default for displayOn = 1.
% The starting vesicle info is taken from entry vindex in the mi.vesicle
% arrays, and the refined information is updated into the mi copy, miNew.
% mode=1 fit amplitude only
% mode=2 fit position and amplitude (default)
% mode=3 fit radius too.
minUnmaskedFraction=0.3; % Don't mask at all if less than this remains of a vesicle.
maxShift=100;  % Å

ndis=size(effCTF,1);
if nargin<9
    displayOn=1;
end;
if nargin<8
    mode=2;  % don't fit radius
end;

mask=single(mask);

if mode>2
    niters=20;  % iterate fitting of the vesicle radius
else
    niters=2;
end;
sigmaT=2;  % SD of prior for shifts
n=size(img,1);
ds=mi.imageSize(1)/n;

% Get the vesicle coordinates relative to ms
vx=(mi.vesicle.x(vindex))/ds+1;  % start with zero-based coordinates
vy=(mi.vesicle.y(vindex))/ds+1;
vr=mi.vesicle.r(vindex,1)/ds;
vs=mi.vesicle.s(vindex,1);

% Get the membrane cross-section density
vd1=mi.vesicleModel;
vd=meDownsampleVesicleModel(vd1,ds)*mi.pixA*ds;

vdOrig=mi0.vesicleModel;
vd0=meDownsampleVesicleModel(vdOrig,ds)*mi0.pixA*ds;

approxLoc=round([vx vy]);
% this is the approximate position of the vesicle relative to
% the center of the full-size image

% Get a centered, cropped portion of the image.
imgc0=ExtractImage(img,approxLoc,ndis);
mask0=ExtractImage(mask,approxLoc,ndis);
% Make a mask for the ccf to prevent excessive shifts
ccmaskR=maxShift/(mi0.pixA*ds)+2;  % radius: add 2 pixels for offset errors
ccmask=fuzzymask(ndis,2,ccmaskR);
% compute the original vesicle for adding back
v=-vs*VesicleFromModel(ndis,vr,vd0,[vx vy]-approxLoc+ndis/2+1);
vv=v(:);
unmaskedFraction=(vv'*(vv.*mask0(:)))/(vv'*vv);
mi.vesicle.af(vindex,1)=unmaskedFraction;  % active fraction

% Do no masking if the majority would be masked.
if unmaskedFraction<minUnmaskedFraction
    mask0=1;
end;
imgc0=imgc0.*mask0;

vfilt=real(ifftn(fftn(v).*effCTF)).*mask0;
% We assume that img has been filtered equivalently to effCTF.  That is, if
% prewhitening has been done to img, it should also be included in effCTF.
imgc=(imgc0+vfilt);  % add back the old vesicle

if displayOn
    subplot(2,2,1);
    imacs(imgc);
    title(vindex(1));
    
    subplot(2,2,2);
    imacs(imgc0);
    title(unmaskedFraction);
    drawnow;
end;

logP=-Radius(ndis).^2./(2*sigmaT^2);  % unscaled log prior
ctr=floor(ndis/2+1);

% initial value for location
F=ones(ndis^2,2);  % we'll use F(:,2) for the constant term.
% const=ones(ndis,ndis);
% F(:,2)=const(:);
P=vr;  % simplex parameter
%
sumT=[ctr ctr];  % accumulate the translations
if mode>2
    vr=Simplex('init',P);
end;
% We update the shift values by locating the CC peak; the amplitude (and dc
% offset) by linear least squares; and if desired, the radius by Simplex.

for iRad=1:niters
    vfit=LLSFit(vr);
    %     v=-VesicleFromModel(ndis,vr,vd,sumT);
    %     fvfilt=fftn(v).*effCTF;
    %     vfilt=real(ifftn(fvfilt)).*mask0;
    %     %     eRef=vfilt(:)'*vfilt(:);  % approx. power in the reference.
    %     % Compute the CCF
    %     cc=fftshift(real(ifftn(fftn(imgc).*conj(fvfilt)))).*ccmask;
    %     logPScaled=logP*hfVar/vs;  % scale up to match ccf
    %     [mxc, xi, yi]=max2di(cc+logPScaled);
    %     sumT=sumT+[xi yi]-ctr;
    %     F(:,1)=vfilt(:);
    %     warning('off','MATLAB:singularMatrix');
    %     a=LinLeastSquares(F,imgc(:));
    %     warning('on','MATLAB:singularMatrix');
    %     vs=a(1);   % vesicle amplitude
    %     vfit=F*a;  % fitted vesicle (amplitude and const)
    
    if mode>2 && iRad<niters
        diff=imgc(:)-vfit;
        err=diff'*diff;
        vr=Simplex(err);  % update the vesicle radius
    end;
end;

if mode>2  % Get the radius from the Simplex centroid, recompute the vesicle,
    %            then do the linear fit one more time.
    vr=Simplex('centroid');
    [vfit,a]=LLSFit(vr);
    %     v=-VesicleFromModel(ndis,vr,vd,sumT);
    %     vfilt=real(ifftn(fftn(v).*effCTF)).*mask0;
    %     cc=fftshift(real(ifftn(fftn(imgc).*conj(fftn(vfilt))))).*ccmask;
    %     logPScaled=logP*hfVar/vs;  % scale up to match ccf
    %     [mxc, xi, yi]=max2di(cc+logPScaled);
    %     sumT=sumT+[xi yi]-ctr;
    %     F(:,1)=vfilt(:);
    %     warning('off','MATLAB:singularMatrix');
    %     a=LinLeastSquares(F,imgc(:));
    %     warning('on','MATLAB:singularMatrix');
    %     vfit=F*a;
end;


diff=imgc(:)-vfit;
vesFit=reshape(vfit,ndis,ndis);
diffIm=reshape(diff,ndis,ndis);

if displayOn
    %     disp(xShift);
    subplot(2,2,4);
    imacs(diffIm);
%     subplot(2,2,3);
%     plot([sect(cc) sect(logPScaled)]);  % plot the CC and the prior
    drawnow;
end;

% update the info structure
miNew=mi;
sh=approxLoc-ndis/2-1;
miNew.vesicle.r(vindex,:)=0;
miNew.vesicle.r(vindex,1)=vr*ds;  % radius isn't changed
miNew.vesicle.x(vindex)=(sumT(1)+sh(1)-1)*ds;  % zero-based position in image
miNew.vesicle.y(vindex)=(sumT(2)+sh(2)-1)*ds;
miNew.vesicle.s(vindex,:)=0;
miNew.vesicle.s(vindex,1)=a(1);


    function [vfit,a]=LLSFit(vr)
        v=-VesicleFromModel(ndis,vr,vd,sumT);
        fvfilt=fftn(v).*effCTF;
        vfilt=real(ifftn(fvfilt)).*mask0;
        %     eRef=vfilt(:)'*vfilt(:);  % approx. power in the reference.
        % Compute the CCF
        cc=fftshift(real(ifftn(fftn(imgc).*conj(fvfilt)))).*ccmask;
        logPScaled=logP*hfVar/vs;  % scale up to match ccf
        [mxc, xi, yi]=max2di(cc+logPScaled);
        sumT=sumT+[xi yi]-ctr;
        F(:,1)=vfilt(:);
        warning('off','MATLAB:singularMatrix');
        a=LinLeastSquares(F,imgc(:));
        warning('on','MATLAB:singularMatrix');
        vfit=F*a;  % fitted vesicle (amplitude and const)
    end


end
