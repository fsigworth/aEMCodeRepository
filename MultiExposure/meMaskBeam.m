% function [img mi]=meMaskBeam(img,mi);

ds=32;  % image downsampling for disc fitting
us1=8;  % padding factor
us2=2;  % padding factor
windowFraction=64; % image width / window width
pa='/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/';
pa='/Users/fred/EMWork/Hideki/150114P/';
% pa='/Users/fred/EMWork/Hideki/140625/Box_AMPAR_Liposome_26_10mMGlu_slot4/';
cd(pa);
miName='Info/001sq02_0_Jan05_12.31.03mi.mat';
miName='Info/sq10_042_Jun25_16.49.5mi.txt';
miName=
% mi=ReadMiFile(miName);
mi=ReadMiFile([pa miName]);
% mi.imageFilenames{1}='sq02_0_001Jan05_12.31.03p4sumAlid1.mrc';
% mi.imageFilenames{1}='sq02_0_001Jan05_12.31.03p4sumAlid2.mrc';

ms=meReadImagesNorm(mi,mi.cpe);
m=ms(:,:,1)+1;
m=ReadEMFile([mi.basePath mi.imagePath mi.imageFilenames{1}]);
n=size(m);
nds=n/ds;
mds=Downsample(m,nds);

figure(1);
SetGrayscale;

imacs(mds);

return
%%
% Fit the overall beam
mb=mds<median(mds(:)/2);  % binary image; 1=dark region
mbp=Crop(mb,nds*us1);
fbp=fftn(mbp);
frame=Crop(SquareWindow(nds,nds(1)/windowFraction),nds*us1);
fFrame=fftn(frame);

% Search for disc
r0=nds(1)/2;
r1=nds(1);
dr=nds(1)/64;
nr=(r1-r0)/dr+1;
cc=zeros([us1*nds,nr],'single');
ccNorm=zeros([us1*nds,nr],'single');
rs=zeros(nr,1);
for ir=1:nr
    r=r0+(ir-1)*dr;
    rs(ir)=r;
    template=1-fuzzymask(nds*us1,2,r,dr);
    fTemplate=conj(fftn(template));
    cc(:,:,ir)=fftshift(real(ifftn(fftn(template).*(fbp))));
    ccNorm(:,:,ir)=fftshift(real(ifftn(fTemplate.*fFrame)));
end;

%%
k=1000;
[val,ix,iy,ir]=max3d(cc./(k+ccNorm));
[ix iy ir]
ctr=nds(1)*us1/2+1;
% ir=25;
% deltaR=-2;
% ix=(ctr-440)*deltaR+ctr;
% iy=(ctr-470)*deltaR+ctr;
% [ix iy ir]


match=1-fuzzymask(nds*us1,2,rs(ir),dr,[ix iy]);
% subplot(221);
% imacs(mbp);

subplot(222);
imacs(single(match)+single(mbp));
title('Overlap of binary image and fit');

subplot(223);
imacs(cc(:,:,ir));

subplot(224);
imacs(ccNorm(:,:,ir));
title('ccNorm');
crossCorr=cc(ix,iy,ir)
normCorr=ccNorm(ix,iy,ir)
subplot(221);
imacs(cc(:,:,ir)./(100+ccNorm(:,:,ir)));
title('Normalized cc');

shifts=[ix iy]-ctr;
r=rs(ir);

return

%%  Now try fitting diffraction pattern
mpad=Crop(mds,nds*us2);   
ctr=nds(1)*us2/2+1;
df=1;
B=10;
alpha=-pi/2;
P0=[shifts r df B alpha];
free=[1 1 1 1 1 0];
p=Simplex('init',P0,free);
frame=Crop(frame,nds*us2);
mpad=mpad.*frame;
niter=300;
py=25;
px1=1;
px2=100;
decayExp=1.5;
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
    subplot(223); imacs(Crop(diff,nds));
    title(num2str(p(3:end)));
    drawnow;
    err=diff(:)'*diff(:);
%%
p=Simplex(err);
end;



