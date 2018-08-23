% function [img mi]=meMaskBeam2(img,mi);
rBeam=2e4;  % angstroms
ds=16;  % image downsampling for disc fitting
us1=8;  % padding factor
us2=2;  % padding factor
windowFraction=64; % image width / window width

% pa='/Users/fred/EMWork/Hideki/140105/Box_AMPAR_Liposome_26_10mMGlu_slot4/';
cd('/Users/fred/EMWork/Hideki/151216/Micrograph')
% cd(pa);
% miName='Info/001sq02_0_Jan05_12.31.03mi.mat';
% mi=ReadMiFile([pa miName]);

% ms=meReadImagesNorm(mi,mi.cpe);
imName='1_016_grid5_Dec16_16.10.52ala.mrc';
imName='1_020_grid5_Dec16_16.13.08ala.mrc';

[m,pixA]=ReadEMFile(imName);
m=RemoveOutliers(m);

n=size(m);
nds=n/ds;
mds=Downsample(m,nds);
% kMask=ones(floor([3838 3710]/16),'single');
% kMask=Crop(kMask,nds);
% masks=zeros([nds 3],'single');
% mask=ones(nds,'single');
% 
% for i=1:3
%     masks(:,:,i)=meGetMask(mi,nds,i);
% end;
figure(1);
SetGrayscale;

imacs(mds);

% return


%%
fhp=.05;
% Fit the overall beam in the downsampled image
% mb=mds<median(mds(:)/2);  % binary image; 1=dark region
mb=mds-mean(mds(:));
mb=GaussHP(mb,fhp);
mbp=Crop(mb,nds*us1);
mbp=GaussHP(mbp,fhp);
fbp=fftn(mbp);
% Create a square padded to us1 times the size
frame=Crop(SquareWindow(nds,nds(1)/windowFraction),nds*us1);
fFrame=fftn(frame);

rs=1.2e4/(pixA*ds);
nr=numel(rs);
dr=nr/20;
cc=zeros([nds*us1 nr],'single');
ccNorm=zeros([nds*us1 nr],'single');

for ir=1:nr
    r=rs(ir);
    template=GaussHP(1-fuzzymask(nds*us1,2,r,dr),fhp);
    fTemplate=conj(fftn(template));
    cc(:,:,ir)=-fftshift(real(ifftn(fftn(template).*(fbp))));
    ccNorm(:,:,ir)=fftshift(real(ifftn(fTemplate.*fFrame)));
end;

% Find the maximum of the normalized cross-correlation
k=1000;
[val,ix,iy,ir]=max3d(cc./(k+ccNorm));
[ix iy ir]
ctr=nds(1)*us1/2+1;
% ir=25;
% deltaR=-2;
% ix=(ctr-440)*deltaR+ctr;
% iy=(ctr-470)*deltaR+ctr;
% [ix iy ir]

drm=1;
match=1-fuzzymask(nds*us1,2,rs(ir),drm,[ix iy]);
% subplot(221);
% imacs(mbp);

subplot(222);
ovl=single(match)+single(mbp);

imacs(ovl);
title('Overlap of binary image and fit');

subplot(223);
imacs(cc(:,:,ir));

subplot(224);
imacs(Crop(ovl,nds));
% imacs(ccNorm(:,:,ir));

title('ccNorm');
crossCorr=cc(ix,iy,ir)
normCorr=ccNorm(ix,iy,ir)
subplot(221);
imacs(cc(:,:,ir)./(k+ccNorm(:,:,ir)));
hold on;
plot(ix,iy,'go');
title('Normalized cc');

shifts=[ix iy]-ctr;
r=rs(ir);

%     subplot(224); imacs(mds);

return

%%  Now try fitting diffraction pattern


% mpad=Crop(mds,nds*us2);
med=median(mds(:));
mdsp=mds;
mdsp(:,1:50)=med;
mdsp(1:50,:)=med;
mdsp=mdsp-mean(mdsp(:));

msk=kMask;
kMask(:,1:50)=0;
kMask(1:50,:)=0;


mpad=Crop(mdsp,nds*us2);

frame=Crop(SquareWindow(nds,nds(1)/windowFraction),nds*us1);

ctr=nds(1)*us2/2+1;
df=1;
B=20;
alpha=-pi/2;
P0=[shifts r df B alpha];
free=[1 1  1 1  0  0];
p=Simplex('init',P0,free);
frame=Crop(frame,nds*us2);
mpad=mpad.*frame;
niter=100;
py=25;
px1=1;
px2=100;
decayExp=1.5;
%     msk=1-BinImage(masks(:,:,3),4);
msk=Crop(kMask,nds*us2);
    for i=1:niter
%%
shifts=p(1:2);
    r=p(3);
    df=p(4);
    B=p(5);
    alpha=p(6);
    disc=fuzzymask(nds*us2,2,r,2,shifts+ctr);
%     disc2=single(disc>.5);
    disc1=single(disc.*frame>.9);
    decay=(RadiusNorm(nds*us2));
%     diffrDisc=real(ifftn(fftn(disc).*ifftshift(exp(-(decay/B).^decayExp).*CTF(nds*us2,1,.025,df,0,0,alpha))));
    diffrDisc=real(ifftn(fftn(disc).*ifftshift(CTF(nds*us2,1,.025,df,0,B,alpha))));
    func=(diffrDisc.*frame.*disc).^2;
    func=(diffrDisc.*frame).^2;
% Subtract model
    amp=func(:)'*mpad(:)/(func(:)'*func(:));
    diff=(mpad-amp*func).*disc1;
    subplot(221); imacs(Crop(mpad,nds));
    subplot(224); plot([mpad(px1:px2,py) amp*func(px1:px2,py)]);

% %     Divide by model
%     m=(mpad./func).*disc1;
%     m(isnan(m))=0;
%     amp=(m(:)'*disc1(:))/(disc1(:)'*disc1(:));
%     diff=m-amp*disc1;
%     subplot(221); imacs(Crop(m,nds));
%     subplot(224); plot([m(px1:px2,py) amp*disc1(px1:px2,py)]);

    subplot(222); imacs(Crop(func,nds));
    subplot(223); imacs(Crop(diff.*msk,nds));
    title(num2str(p(3:end)));
    drawnow;
    err=diff(:)'*diff(:);
%%
p=Simplex(err);
end;

%% Show 1D display
mdsz=mds;
mdsz(:,1:5)=0;

mx=Crop(mdsz,size(mdsz)*1.5);
subplot(222);
imacs(mx)

mxr=grotate(mx,-.85);
subplot(221);
imacs(mxr);

mx1=mean(mxr(:,170:200),2);
subplot(223);
plot(mx1)

fux=Crop((amp)*func,size(mdsz)*1.5);
fxr=grotate(fux,-.85);
subplot(224);
imacs(fxr);

fx1=mean(fxr(:,170:200),2);

mx1=mean(mxr(:,170:200),2);
subplot(223);
plot([mx1 fx1])

