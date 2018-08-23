function rsQuickFitVesicleGeneral
z=load('env');
effCTF=z.effCTF;
ds=z.ds;
m=z.m;
mi=z.mi;
mi1=z.mi1;
ind=z.ind;
ndis=z.ndis;
pixA=z.pixA;
mmask=z.mmask;
ms=z.ms;
msub=z.msub;
msubf=z.msubf;
miRef=load('/Users/fred/EMWork/Hideki/140411/Kvbeta_frB2_liposome_offsetneg2.5/Info/sq03_1_001mi.mat');
miRef=miRef.mi;
img=msubf;
mask=mmask;
mi=mi1;
mi0=mi1;
vIndex=ind;
mode=3;
doDisplay=1;
a0=-0.0;
mi.ctf(1).alpha=a0;
mi.ctf(2).alpha=a0;
sigmaN=.001;
effCTF=ifftshift(meGetEffectiveCTF(mi,ndis,ds));

figure(1);
SetGrayscale;

% function [miNew, diffIm, vesFit]=rsQuickFitVesicleDeformed(img,mask,mi,mi0,vindex,effCTF,hfVar,mode,displayOn)
% function [miNew,diffIm,vesFit]=rsQuickFitVesicleDeformed(img,mask,mi,mi0,vindex,effCTF,hfVar,mode,displayOn)
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
deltaR=1;
nTerms=4;

n2=ndis^2;

ndis=size(effCTF,1);

mask=single(mask);

n=size(img,1);
ds=mi.imageSize(1)/n;

% Get the vesicle coordinates relative to ms
vx=(mi.vesicle.x(vIndex))/ds+1;  % start with zero-based coordinates
vy=(mi.vesicle.y(vIndex))/ds+1;
vr=mi.vesicle.r(vIndex)/ds;
vs=mi.vesicle.s(vIndex);


% Get the membrane cross-section density
% Density for adding a vesicle back
vdOrig=mi0.vesicleModel;
vd0=meDownsampleVesicleModel(vdOrig,ds)*mi0.pixA*ds;

% % Density from the new mi
% vd1=mi.vesicleModel;
% vd=meDownsampleVesicleModel(vd1,ds)*mi.pixA*ds;

% Use the ref vesicle model instead
dsm=pixA/miRef.pixA;
vd=meDownsampleVesicleModel(miRef.vesicleModel,dsm)*pixA;

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
org=[vx vy]-approxLoc+ndis/2+1; % Origin in extracted image
v0=VesicleFromModel(ndis,vr,vd0,org);

% vFit=vs*v0;
% 
% 
% v00=VesicleFromModelDeformed(ndis,vr,vd0,org);
% Somehow, VesicleFromModelDeformed(ndis,vr-1) gives the same result as
% VesicleFromModel(ndis,vr)???
v=-vs*v0;
vv=v(:);
unmaskedFraction=(vv'*(vv.*mask0(:)))/(vv'*vv);
% Do no masking if the majority would be masked.
if unmaskedFraction<minUnmaskedFraction
    mask0=1;
end;
% nIters=10;
% for j=1:nIters
% make the derivative wrt radius
deltas=zeros(1,nTerms);
deltas(1)=deltaR;
dv0=(-v0+VesicleFromModelGeneral(ndis,vr+deltas,vd,org))/(deltaR);

fvFilt=fftn(v0).*effCTF;
vFilt=real(ifftn(fvFilt)).*mask0;
fvdFilt=fftn(dv0).*effCTF;
vdFilt=real(ifftn(fvdFilt)).*mask0;


%%%-------------------------test code
% imgu=-vs*VesicleFromModelDeformed(ndis,[vr .5i .2 0],vd0,org,[1 .0 .0]);
% % est. a in this case is [-.001(=scaler) ..... -1e-5(=scaler * value/3) x 3]
% imgc0=real(ifftn(fftn(imgu).*effCTF))+sigmaN*randn(ndis);

% imgc=(imgc0-vFilt).*mask0;  % masked, subtracted image


% [r2D,theta2D]=Radius(ndis,org);
% 
% theta=reshape(theta2D,n2,1);
% w=exp(1i*theta);  % complex exponential
% wi=ones([n2 nTerms],'single');
% % Get powers of the complex exponential
% for j=2:nTerms
%     wi(:,j)=wi(:,j-1).*w(:);
% end;

vFitU=0;

vMask=fuzzymask(ndis,2,vr+100/pixA,50/pixA,org)...
     -fuzzymask(ndis,2,vr-60/pixA,50/pixA,org);
imgc1=imgc0.*vMask;

subplot(2,2,1);
imacs(imgc0);

x0=zeros(1,4*nTerms-2);
x0(1)=vr;
x0(2*nTerms-1)=1;

opts=optimoptions('lsqnonlin','TolFun',1e-6,'TolX',1e-5,...
    'FinDiffRelStep',.01,'DiffMinChange',1e-5,'Display','off');
[a,resNorm,residual,flag,output]=lsqnonlin(@FitFcn,x0,[],[],opts);
output.funcCount
subplot(2,2,3);
imacs(imgc0-vFitU);

[rs,ws]=InsertParameters(a);

% Make sure the mi structure has the requisite fields
if ~isfield(mi.vesicle,'rx')
    mi.vesicle.rx=mi.vesicle.r;
    mi.vesicle.wx=ones(size(mi.vesicle.s),'single');
end;
nrs=numel(rs);
if size(mi.vesicle.rx,2)<nrs
    mi.vesicle.rx(1,nrs)=0;  % expand the array
    mi.vesicle.wx(1,nrs)=0;
end;

mi.vesicle.rx(vIndex,:)=rs*ds;
mi.vesicle.wx(vIndex,:)=ws;

return

% function for lsqnonlin()
    function F=FitFcn(x)
%         disp(x);
        [rs,ws]=InsertParameters(x);
        vFit0=vs*VesicleFromModelGeneral(ndis,rs,vd,org,ws);
        vFitU=real(ifftn(fftn(vFit0).*effCTF));
        vFit=vFitU.*vMask;
        F=double(imgc1-vFit);
        subplot(2,2,2);
        imacs(F);
        subplot(2,2,4);
        imacs(vFit);
        subplot(2,2,3);
        imacs(imgc0-vFitU);
        drawnow;
    end


    function [rs,ws]=InsertParameters(x)
        nt=(numel(x)+2)/4;
        rs=[x(1) .01*x(2:nt) + .01i*x(nt+1:2*nt-1)];
        ws=[x(2*nt) .01*x(2*nt+1:3*nt-1) + .01i*x(3*nt:4*nt-2)];
    end




end


% 
% 
% % Set up the LS problem
% F=ones(ndis^2,(1+deltaRadius)*nTerms+1)*100;
% for j=1:nTerms
%     F(:,j)=vFilt(:).*wi(:,j);  % amplitudes
%     if deltaRadius
%         F(:,j+nTerms)=vdFilt(:).*wi(:,j); % derivatives
%     end;
% end;
% % last term is constant
% imgc=imgc0;
% % imgc=vFilt+randn(ndis)*10;
% a=LinLeastSquares(F,double(imgc(:)));
% 
% vr(2:nTerms)=a(nTerms+2:2*nTerms)/a(1)
% %%
% afit=real(F*a);
% afit=reshape(afit,ndis,ndis);
% 
% subplot(2,2,1);
% imacs(imgc);
% 
% subplot(2,2,4);
% imacs(afit);
% 
% subplot(2,2,2);
% imacs(imgc-afit);
% 
% subplot(2,2,3);
% imacs(GaussFilt(imgc-afit,.2));
% 
% a
% return





% 
% 
% vfilt=real(ifftn(fftn(v).*effCTF)).*mask0;
% % We assume that img has been filtered equivalently to effCTF.  That is, if
% % prewhitening has been done to img, it should also be included in effCTF.
% imgc=(imgc0+vfilt);  % add back the old vesicle
% 
% 
% 
% 
% 
% if displayOn
%     subplot(2,2,1);
%     imacs(imgc);
%     title(vindex(1));
%     
%     subplot(2,2,2);
%     imacs(imgc0);
%     title(unmaskedFraction);
%     drawnow;
% end;
% 
% % logP=-Radius(ndis).^2./(2*sigmaT^2);  % unscaled log prior
% % ctr=floor(ndis/2+1);
% %
% % % initial value for location
% % F=ones(ndis^2,2);  % we'll use F(:,2) for the constant term.
% % % const=ones(ndis,ndis);
% % % F(:,2)=const(:);
% % P=vr;  % simplex parameter
% % %
% % sumT=[ctr ctr];  % accumulate the translations
% % if mode>2
% %     vr=Simplex('init',P);
% % end;
% % % We update the shift values by locating the CC peak; the amplitude (and dc
% % % offset) by linear least squares; and if desired, the radius by Simplex.
% %
% % for iRad=1:niters
% %     vfit=LLSFit(vr);
% %     %     v=-VesicleFromModel(ndis,vr,vd,sumT);
% %     %     fvfilt=fftn(v).*effCTF;
% %     %     vfilt=real(ifftn(fvfilt)).*mask0;
% %     %     %     eRef=vfilt(:)'*vfilt(:);  % approx. power in the reference.
% %     %     % Compute the CCF
% %     %     cc=fftshift(real(ifftn(fftn(imgc).*conj(fvfilt)))).*ccmask;
% %     %     logPScaled=logP*hfVar/vs;  % scale up to match ccf
% %     %     [mxc, xi, yi]=max2di(cc+logPScaled);
% %     %     sumT=sumT+[xi yi]-ctr;
% %     %     F(:,1)=vfilt(:);
% %     %     warning('off','MATLAB:singularMatrix');
% %     %     a=LinLeastSquares(F,imgc(:));
% %     %     warning('on','MATLAB:singularMatrix');
% %     %     vs=a(1);   % vesicle amplitude
% %     %     vfit=F*a;  % fitted vesicle (amplitude and const)
% %
% %     if mode>2 && iRad<niters
% %         diff=imgc(:)-vfit;
% %         err=diff'*diff;
% %         vr=Simplex(err);  % update the vesicle radius
% %     end;
% % end;
% %
% % if mode>2  % Get the radius from the Simplex centroid, recompute the vesicle,
% %     %            then do the linear fit one more time.
% %     vr=Simplex('centroid');
% %     [vfit,a]=LLSFit(vr);
% %     %     v=-VesicleFromModel(ndis,vr,vd,sumT);
% %     %     vfilt=real(ifftn(fftn(v).*effCTF)).*mask0;
% %     %     cc=fftshift(real(ifftn(fftn(imgc).*conj(fftn(vfilt))))).*ccmask;
% %     %     logPScaled=logP*hfVar/vs;  % scale up to match ccf
% %     %     [mxc, xi, yi]=max2di(cc+logPScaled);
% %     %     sumT=sumT+[xi yi]-ctr;
% %     %     F(:,1)=vfilt(:);
% %     %     warning('off','MATLAB:singularMatrix');
% %     %     a=LinLeastSquares(F,imgc(:));
% %     %     warning('on','MATLAB:singularMatrix');
% %     %     vfit=F*a;
% % end;
% %
% %
% % diff=imgc(:)-vfit;
% % vesFit=reshape(vfit,ndis,ndis);
% % diffIm=reshape(diff,ndis,ndis);
% %
% % if displayOn
% %     %     disp(xShift);
% %     subplot(2,2,4);
% %     imacs(diffIm);
% %     subplot(2,2,3);
% %     %     imacs(vesFit);
% %     plot([sect(cc) sect(logPScaled)]);  % plot the CC and the prior
% %     drawnow;
% % end;
% %
% % % update the info structure
% % miNew=mi;
% % sh=approxLoc-ndis/2-1;
% % miNew.vesicle.r(vindex)=vr*ds;  % radius isn't changed
% miNew.vesicle.x(vindex)=(sumT(1)+sh(1)-1)*ds;  % zero-based position in image
% miNew.vesicle.y(vindex)=(sumT(2)+sh(2)-1)*ds;
% miNew.vesicle.s(vindex)=a(1);
%
%
%     function [vfit,a]=LLSFit(vr)
%         v=-VesicleFromModel(ndis,vr,vd,sumT);
%         fvfilt=fftn(v).*effCTF;
%         vfilt=real(ifftn(fvfilt)).*mask0;
%         %     eRef=vfilt(:)'*vfilt(:);  % approx. power in the reference.
%         % Compute the CCF
%         cc=fftshift(real(ifftn(fftn(imgc).*conj(fvfilt)))).*ccmask;
%         logPScaled=logP*hfVar/vs;  % scale up to match ccf
%         [mxc, xi, yi]=max2di(cc+logPScaled);
%         sumT=sumT+[xi yi]-ctr;
%         f=ones(ndis^2,2);
%         f(:,1)=vfilt(:);
%         warning('off','MATLAB:singularMatrix');
%         a=LinLeastSquares(F,imgc(:));
%         warning('on','MATLAB:singularMatrix');
%         vfit=f*a;  % fitted vesicle (amplitude and const)
%     end
%
%
% end
