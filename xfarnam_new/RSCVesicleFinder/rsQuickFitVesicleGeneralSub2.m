function [miNew, fitIm, vesFit]=rsQuickFitVesicleGeneralSub2(img,mask,mi,mi0,vIndex,effCTF,hfVar,mode,nTerms,displayOn)
% function [miNew, fitIm, vesFit]=rsQuickFitVesicleGeneralSub2(img,mask,mi,mi0,vIndex,effCTF,hfVar,mode,nTerms,displayOn)
% Version 2 uses an initial Simplex fit to get close.
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

nCallsMin=100;  % number of fcn calls before starting to show display

% Rescaling factors for fit
sc0=.01;   % for higher-order geometry factors
scw=.001;  % for amplitude factor

% Check to see if any higher-order terms will be used.
if nTerms<3
    [miNew, fitIm, vesFit]=rsQuickFitVesicle(img,mask,mi,mi0,vIndex,effCTF,hfVar,mode,displayOn);
    return
end;

minUnmaskedFraction=0.3; % Don't mask at all if less than this remains of a vesicle.

ndis=size(effCTF,1);    % size of image we'll display (and fit)
cdis=ceil((ndis+1)/2);  % center coordinate

mask=single(mask);

n=size(img,1);          % Downsampled micrograph
ds=mi.imageSize(1)/n;
pixA=mi.pixA*ds;

% Get the vesicle coordinates relative to img
vx=(mi.vesicle.x(vIndex))/ds+1;  % start with zero-based coordinates
vy=(mi.vesicle.y(vIndex))/ds+1;
vr=mi.vesicle.r(vIndex,1)/ds;
vs=mi.vesicle.s(vIndex,1);


% Get the membrane cross-section density
% -the density for adding a vesicle back from past subtraction
vdOrig=mi0.vesicleModel;
vd0=meDownsampleVesicleModel(vdOrig,ds)*mi0.pixA*ds;

% Density from the new mi in case it's different
vd1Orig=mi.vesicleModel;
vd1=meDownsampleVesicleModel(vd1Orig,ds)*mi.pixA*ds;

approxLoc=round([vx vy]);
% this is the approximate position of the vesicle relative to
% the center of the full-size image

% Get a centered, cropped portion of the image.
imgC0=ExtractImage(img,approxLoc,ndis);
org0=[vx vy]-approxLoc; % relative origin in extracted image; ndis must be even
v0=VesicleFromModelGeneral(ndis,vr,vd0,org0+cdis,vs); % existing vesicle model
vFilt0=-real(ifftn(fftn(v0).*effCTF));  % filtered old vesicle to add back

% Get the extracted mask
mask0=ExtractImage(mask,approxLoc,ndis);
vv=v0(:);  % compute the overlap of raw vesicle model with the mask.
unmaskedFraction=(vv'*(vv.*mask0(:)))/(vv'*vv);
mi.vesicle.af(vIndex,1)=unmaskedFraction;
% Do no masking if the majority would be masked.
if unmaskedFraction<minUnmaskedFraction
    mask0=1;
end;

% Fit only over an annular region
vMask=fuzzymask(ndis,2,vr+80/pixA,50/pixA,org0+cdis)...
    -fuzzymask(ndis,2,vr-60/pixA,50/pixA,org0+cdis);
imgR1=(imgC0+vFilt0).*mask0;  % Basic mask, vesicle-restored
imgToFit=imgR1.*vMask; % Restored, masked image

if displayOn
    subplot(2,2,1);
    imacs(imgR1);  % display the restored image without annular mask.
    title(vIndex);
    subplot(2,2,3);
    imacs(imgC0);
end;

% Do the fitting.

% Set the initial parameter values a0, inverse of DecodeParameters(a)
%  a0(1:2nTerms-1), total of 2nTerms-1 values
%    real r values 0, 2...nTerms; imag r values 2...nTerms
%  a0(2nTerms:4nTerms), total of 2nTerms+1 values
%    real weights 0...nTerms; imag weights 1...nterms
%  a0(4nTerms+1:4nTerms+2):  orgX, orgY.
% for a total of (nTerms+nTerms-1)+(nTerms+1+nTerms) + 2 = 4nTerms+2
% **Except in the case nTerms=1, then a has 2+2 values.
if nTerms<2 % case of nTerms=1; create a 1x4 vector
    a0=[vr vs/scw org0];
    msk=logical([1 0 1 1]);
else
    a0=zeros(1,4*nTerms+2);
    a0(1)=vr;  % zero-order radius
    a0(2*nTerms)=vs/scw; % zero-order weight
    a0(4*nTerms+1:4*nTerms+2)=org0; % shift x and y
    msk=false(1,4*nTerms+2);
    msk(1:2*nT1-1)=true;  % fit only the first nT1 radius terms
    msk(4*nTerms+1:4*nterms+2)=true;    
end;

% Get the complex exponentials
[R,theta]=Radius2(ndis,org0+cdis);
n2=prod(ndis);
w=exp(1i*theta(:));  % complex exponential
wi=ones([n2 nterms],'single');
% Get powers of the complex exponential
for j=2:nterms
    wi(:,j)=wi(:,j-1).*w(:);
end;

vesFit=zeros(ndis,ndis,'single');  % place to put the unmasked fitted vesicle image.

% Simplex fitting is done here-------

msk=false(size(a0));
a1=Simplex('init',a0,steps,msk);
for si=1:nsi
    [rs,ws,org]=DecodeParameters(a)
    v0=-VesicleFromModelGeneral(ndis,rs,vd,org);
    F=




opts=optimoptions('lsqnonlin','TolFun',1e-6,'TolX',1e-4,...
    'FinDiffRelStep',.01,'DiffMinChange',1e-5,'Display','off');
% [a,resNorm,residual,flag,output]=lsqnonlin(@FitFcn,a0,[],[],opts);
fitOk=1;
try
a=lsqnonlin(@FitFcn,a0,[],[],opts);
catch
    a=a0;
    fitOk=0;
end;
if displayOn
    subplot(2,2,3);
    imacs((imgR1-vesFit));  % unmasked final fit
end;

% Get the returned values
FitFcn(a);  % This sets the global vesFit, which is returned
fitIm=imgR1;

% Make sure the mi structure has the requisite fields.  Note that we
% asssume that vesicle.r and vesicle.s are the same size always.
miNew=mi;
szr=size(miNew.vesicle.r,2);
if szr<nTerms+1
    miNew.vesicle.r(1,nTerms+1)=0;  % expand the arrays
    miNew.vesicle.s(1,nTerms+1)=0;
elseif szr>nTerms+1
    miNew.vesicle.r(vIndex,:)=0;  % Make sure no left-over elements
    miNew.vesicle.s(vIndex,:)=0;
end;
[rs,ws,org]=DecodeParameters(a);
ws=ws*fitOk;  % NaNs if it's a bad fit.
miNew.vesicle.r(vIndex,1:nTerms+1)=rs*ds;  % now a vector
miNew.vesicle.s(vIndex,1:nTerms+1)=ws;
miNew.vesicle.x(vIndex)=(org(1)+approxLoc(1)-1)*ds;
miNew.vesicle.y(vIndex)=(org(2)+approxLoc(2)-1)*ds;

return

% function for lsqnonlin()
    function F=FitFcn(a)
        %         Returns img-fit given the parameter vector a.
        
        nCalls=nCalls+1;
        [rs,ws,org]=DecodeParameters(a);
        vFit0=VesicleFromModelGeneral(ndis,rs,vd1,org+cdis,ws);
        vesFit=-real(ifftn(fftn(vFit0).*effCTF)).*mask0;
        vFit=vesFit.*vMask;
        F=double(imgToFit-vFit);
        if displayOn && nCalls>nCallsMin && mod(nCalls,10)==0
            %             subplot(2,2,2);
            %             imacs(F);
            subplot(2,2,4);
            imacs(vFit);
            subplot(2,2,3);
            imacs(imgR1-vesFit);
            title(nCalls);
            drawnow;
        end;
    end

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


    function [rs,ws,org]=DecodeParameters(a)
        nt=(numel(a)-2)/4;
        if nt<1
            rs=a(1);
            ws=a(2);
            org=a(3:4);
        else
            %             rs and ws will be complex vectors 1 x (nt+1)
            %        Zero-order terms rs(1) and ws(1) are real.
            %        The radius expansion skips the 1st order term (2nd vector element)
            rs=[a(1) 0 (sc0*a(2:nt) + sc0*1i*a(nt+1:2*nt-1))];
            %         Weight expansion does not skip the 1st order term.
            ws=[a(2*nt)*scw sc0*a(2*nt+1:3*nt) + sc0*1i*a(3*nt+1:4*nt)];
            org=a(4*nt+1:4*nt+2);
        end;
    end




end
