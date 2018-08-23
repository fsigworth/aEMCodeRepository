function miNew=rsFineFitVesicle2(img,mi,mi0,vindex,ndis,displayOn,img0)
% version 2 uses the built-in Levinberg-Marquardt algorithm.
% function miNew=rsFineFitVesicle2(img,mi,vindex,ndis,displayOn)
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
showOriginal = nargin >= 6;

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
nd=ceil((numel(vd1)-1)/ds+1);
vd=Interpolate1(vd1,nd,1/ds)*ds*mi.pixA; % scale by overall pixel size.

vdOrig=mi0.vesicleModel;
nd=ceil((numel(vdOrig)-1)/ds+1);
vd0=Interpolate1(vdOrig,nd,1/ds)*ds*mi0.pixA; % scale by overall pixel size.


approxShift=round([vx vy]-n/2);
% this is the approximate position of the vesicle relative to
% the center of the full-size image

% Get a centered, cropped portion of the image.
imgc=Crop(circshift(img,-approxShift),ndis);
if displayOn
    subplot(2,2,2);
    imacs(imgc);
    drawnow;
end;

% 
% 
% 
% 
% 
% if showOriginal && displayOn
%     img0c=Crop(circshift(img0,-approxShift),ndis);
%     subplot(2,2,1);
%     imacs(img0c);
%     drawnow;
% end;

% Get the filter function
H=ifftshift(meGetEffectiveCTF(mi,ndis,ds));

% add back the original vesicle
    v=-vs*VesicleFromModel(ndis,vr,vd0,[vx vy]-approxShift+(ndis-n)/2);
    vfilt=real(ifftn(fftn(v).*H));
    imgc=imgc+vfilt;

    if displayOn
    subplot(2,2,1);
    imacs(imgc);
    drawnow;
end;


P0=double([vr ndis/2+1 ndis/2+1]);

iter=0;
        options=optimset('TolX',1e-4,'FinDiffRelStep',1e-5,'Display','off');
%             'Diagnostics','off');

[P,resnorm,residual,exitFlag]=lsqnonlin(@GetResidual,P0,0,inf,options);

% update the info structure
miNew=mi;
sh=approxShift-ndis/2+n/2;
miNew.vesicle.r(vindex)=P(1)*ds;
miNew.vesicle.x(vindex)=(P(2)+sh(1)-1)*ds;  % zero-based position in image
miNew.vesicle.y(vindex)=(P(3)+sh(2)-1)*ds;
miNew.vesicle.s(vindex)=a(1);

% exitFlag

% disp([P(1) P(2) P(3) a']);  % show the main fit


    function F=GetResidual(x)
        iter=iter+1;
    r=x(1);
    v=-VesicleFromModel(ndis,r,vd,x(2:3));
    vfilt=real(ifftn(fftn(v).*H));
    %     vpow=vfilt(:)'*vfilt(:);
    const=ones(ndis,ndis)*100;
    R=[vfilt(:) const(:)];
    a=LinLeastSquares(R,imgc(:));
    vfit=R*a;
    F=double(imgc(:)-vfit);
    
    if mod(iter,25)==0 && displayOn
        subplot(2,2,4);
        imacs(reshape(F,ndis,ndis));
        subplot(2,2,3);
        imacs(vfilt);
        title(iter);
        drawnow;
    end;
%     
    end
end
