% RotEMTRPV1.m
% E-M algorithm for restoring rotated, noisy images.


rotAngle=45; % phi angle for TRPV1 tetramer
rotAngle=0
cropFactor=2;
psiAngle=180;
cd ~/aEMCodeRepository/EMClass/
[m,s]=ReadMRC('emd_5778.map');
if mod(rotAngle,90)~=0
    disp('Rotating...');
    mr=rsRotateImage(m,rotAngle); % rotate about z axis
    disp('done.');
else
    mr=m;
end;
mr=rsRotateImage(squeeze(sum(mr,1)),psiAngle); % projection
m0=Crop(mr,512);
%%
m12=Downsample(m0,256);
pixA=s.pixA*2;

n=size(m12,1);
ndis=3*n/8;
figure(4);

kV=300;
def0=1;
def1=2;

Cs=2;
B=100/4;
% B=0;
alpha=.05;
deltadef=0;
theta=pi/4;

m1=Crop(m12,ndis);
mc=m1/std(m1(:));

nDefs=20;
nPerDef=5000;

defStep=(def1-def0)/nDefs;
defs=def0:defStep:def1-defStep;
% Make dataset
sigma=20;
rImgAng=120;
showImgs=1;

dImgs=zeros(ndis,ndis,nPerDef,nDefs,'single');
cts=zeros(ndis,ndis,nDefs,'single');
A0=rsRotateImage(mc,rImgAng);
for iDef=1:nDefs
    cts(:,:,iDef)=ifftshift(-CTF(ndis,pixA,EWavelength(kV),defs(iDef),Cs,B,alpha,deltadef,theta));
    rImg0=real(ifftn(fftn(A0).*cts(:,:,iDef)));
    if showImgs
        fc=.2;
        subplot(4,4,3);
        imags(GaussFilt(rImg0,fc));
        axis off;
        subplot(4,4,4);
        title(defs(iDef));
        imags(GaussFilt(rImg0+randn(ndis)*sigma,fc));
        axis off;
        drawnow;
    end;
    
    for iPer=1:nPerDef
        N1=(randn(ndis)*sigma);
        dImgs(:,:,iPer,iDef)=rImg0+N1;
    end;
end;

%% ---EM iterations--   
    iter=0;
    nIters=10;
ref0=dImgs(:,:,1);
ref0=mc+randn(ndis)*sigma;


% rAngs=0:2:360;
rAngs=100:2:140;
nRefs=numel(rAngs);
rRefDs=zeros(ndis,ndis,nRefs,'single');
A=zeros(ndis,ndis);
k=1e5;
%%
while iter<nIters
    iter=iter+1
    rRefs=rsRotateImage(ref0,rAngs);
    cfDefAccum=zeros(ndis,ndis);
    defAccum=zeros(ndis,ndis);
    c2Accum=zeros(ndis,ndis);
        SS=zeros(nRefs,1);
    
     for iDef=1:nDefs % Loop over defocus
%       for iDef=1:3
        for iRef=1:nRefs
        rRefDs(:,:,iRef)=real(ifftn(fftn(rRefs(:,:,iRef)).*cts(:,:,iDef)));
        end;
        dAccums=zeros(ndis,ndis,nRefs);
        for iPer=1:nPerDef % loop over group
            img=dImgs(:,:,iPer,iDef);
            for iRef=1:nRefs % over refs
                dif=img-rRefDs(:,:,iRef);
                SS(iRef)=-(dif(:)'*dif(:))/(2*sigma);
            end;
            mxSS=max(SS);
            SS=SS-mxSS;
            pSS=exp(SS);
            norm=sum(pSS);
            pSS=pSS/norm; % P(psi);
            for iRef=1:nRefs
                dAccums(:,:,iRef)=dAccums(:,:,iRef)+pSS(iRef)*img;
            end;
        end;
        % now we can de-rotate the accumulators
        defAccum=sum(rsRotateImage(dAccums,-rAngs),3);
        cfDefAccum=cfDefAccum+fftn(defAccum).*cts(iDef);
        c2Accum=c2Accum+nPerDef*cts(:,:,iDef).^2;

        
        subplot(4,4,4);
imags(defAccum);
        subplot(121);
imagsarray(dAccums);
title(iDef,'color','w');
drawnow;
      
    end;
    %%
    A1=real(ifftn(cfDefAccum./(c2Accum+k)));
    subplot(443);
    imags(A1);
    A(:,:,iter)=A1;
    ref0=A1;
    
end;

    
