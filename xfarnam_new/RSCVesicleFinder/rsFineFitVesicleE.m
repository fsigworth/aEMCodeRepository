function [miNew diffIm vesFit]=rsFineFitVesicleE(img,mi,mi0,vindex,ndis,displayOn,mx)
% function miNew=meFineFitVesicle(img,mi,vindex,ndis,displayOn)
% function [vesfit info]=FitOneVesicle(sirv,vesinfo,ginfo,Feff)
% Given a subtracted image with one vesicle restored, fit by
% least-squares a vesicle model with fractional shifts and radius.  The
% starting vesicle info is taken from entry vindex in the mi.vesicle
% arrays, and the refined information is written back into this place.
% A small image of size ndis is extracted from img and this is the portion
% that is fitted.  The default value for ndis = 128 and the default for
% displayOn = 1.

if nargin<4
    ndis=128;
end;
if nargin<5
    displayOn=1;
end;

niters=401;
n=size(img,1);
ds=mi.imageSize(1)/n;
% pixA=mi.pixA*ds;
% Get the coordinates and radius, scaled down by ds
vx=(mi.vesicle.x(vindex))/ds+1;  % assume zero-based coordinates
vy=(mi.vesicle.y(vindex))/ds+1;
vr=mi.vesicle.r(vindex)/ds;
vs=mi.vesicle.s(vindex);

% Get the membrane cross-section density
% Get the membrane cross-section density
vd1=mi.vesicleModel;
% nd=ceil((numel(vd1)-1)/ds+1);
% vd=Interpolate1(vd1,nd,1/ds)*ds*mi.pixA; % scale by overall pixel size.
vd=meDownsampleVesicleModel(vd1,ds)*mi.pixA*ds;

vdOrig=mi0.vesicleModel;
vd0=meDownsampleVesicleModel(vdOrig,ds)*mi0.pixA*ds;
% nd=ceil((numel(vdOrig)-1)/ds+1);
% vd0=Interpolate1(vdOrig,nd,1/ds)*ds*mi0.pixA; % scale by overall pixel size.

approxShift=round([vx vy]-n/2);
% this is the approximate position of the vesicle relative to
% the center of the full-size image

% Get a centered, cropped portion of the image.
imgc0=Crop(circshift(img,-approxShift),ndis);
if displayOn
    subplot(2,2,2);
    imacs(imgc0);
    drawnow;
end;

% Get the filter function
H=ifftshift(meGetEffectiveCTF(mi,ndis,ds));

% add back the original vesicle
     v=-vs*VesicleFromModel(ndis,vr,vd0,[vx vy]-approxShift+(ndis-n)/2);
     vfilt=real(ifftn(fftn(v).*H));
    imgc=imgc0+vfilt;

    if displayOn
    subplot(2,2,1);
    imacs(imgc);
    title(vindex(1));
    drawnow;
end;


stepsizes=[.2 .2 .2 .01 .01 .01 .01];
P=[vr ndis/2+1 ndis/2+1 0 0 0 0];

P=Simplex('init',P,stepsizes);

for iter=1:niters
    % Construct the vesicle and interior functions
    r=P(1);
    e=P(4:7);
    v=-EllipticalVesicleFromModel(ndis,r,e,vd,P(2:3));
    vfilt=real(ifftn(fftn(v).*H));
    %     vpow=vfilt(:)'*vfilt(:);
    const=ones(ndis,ndis)*100;
    F=[vfilt(:) const(:)];
    a=LinLeastSquares(F,imgc(:));
    vfit=F*a;
    diff=imgc(:)-vfit;
    err=diff'*diff;
    
    vesFit=reshape(vfit,ndis,ndis);
    diffIm=reshape(diff,ndis,ndis);
    %     subplot(2,2,3);
    %     imacs(vesfit);
    
    if mod(iter,50)==1 && displayOn
        subplot(2,2,4);
        imacs(diffIm);
        subplot(2,2,3);
        imacs(vesFit);
        title(iter);
        drawnow;
    end;
    
    P=Simplex(err);
end;

P=Simplex('centroid');
%% Construct the model image
    r=P(1);
        e=P(4:7);

    v=-EllipticalVesicleFromModel(ndis,r,e,vd,P(2:3));
    vfilt=real(ifftn(fftn(v).*H));
    F=[vfilt(:) const(:)];
    a=LinLeastSquares(F,imgc(:));
    vfit=F*a;  % Ignore the constant term
    diff=imgc(:)-vfit;
    vesFit=reshape(vfit,ndis,ndis);
    diffIm=reshape(diff,ndis,ndis);  % residual image
    
% update the info structure
miNew=mi;
sh=approxShift-ndis/2+n/2;
miNew.vesicle.r(vindex)=P(1)*ds;
miNew.vesicle.x(vindex)=(P(2)+sh(1)-1)*ds;  % zero-based position in image
miNew.vesicle.y(vindex)=(P(3)+sh(2)-1)*ds;
miNew.vesicle.s(vindex)=a(1);

if nargout>1  % We do the angular ananalysis
    %%%%%%%%%%%%%%
    
    
end;




% disp([P(1) P(2) P(3) a']);  % show the main fit

end
