% meEvalRun
ds=2;
i1=1;
i2=2;
% original sigma^2 is 5 (10e exposures) or 10 (20e exposure).
noiseSigma=0;
n=size(m,1);
m1=AffineTransform(BinImage(m(:,:,i1)+noiseSigma*randn(n,n),ds),(mi.mergeMatrix(:,:,i1)));
m2=AffineTransform(BinImage(m(:,:,i2)+noiseSigma*randn(n,n),ds),(mi.mergeMatrix(:,:,i2)));
% m1=BinImage(m(:,:,i1),2);
% m2=BinImage(m(:,:,i2),2);

np=size(m1,1);
pixA=mi.pixA*ds;
cu1=single(CTF(np,pixA,mi.ctf(i1)));
cu2=single(CTF(np,pixA,mi.ctf(i2)));
radexp=0;
pw=Radius(np).^radexp;
pw(np/2+1,np/2+1)=0;

filt=fftshift(sign(cu1).*sign(cu2).*(abs(cu1.*cu2)).*pw);


% m1=BinImage(m(:,:,1),2);
% m2=BinImage(m(:,:,3),2);

pixA=mi.pixA*ds;
ndis=64;
fig=5;
CPars=mi.ctf(i1);
CPars(2)=mi.ctf(i2);
npanels=8;

cc=fftshift(real(ifftn(fftn(m1).*conj(fftn(m2)).*filt)));
subplot(1,3,1);
imacs(Crop(cc,128));

subplot(132);
plot(sect(Crop(cc,1024)));
[val i j]=max2di(Crop(cc,128))

drawnow;

% return

%%

% meEvalLocalCorrelations2(m1,m2,pixA,CPars,npanels,ndis,fig);

% function [dxs dys ccs]=meEvalLocalCorrelations2(m1,m2c,pixA,CPars,npanels,ndis,fig)
% function [dxs dys ccs]=meEvalLocalCorrelations(m1,m2c,CPars,npanels,ndis,fig)
% For two aligned images for merging, compute local cross-correlations
% and return the displacements in each of npanels x npanels parts of the
% images.  CPars is a two-element array of the CTF parameter structures.
% ndis the the number of pixels to return in the cross-correlation images
% ccs.
% This script can be run after meMergeExposures to check the quality of
% local alignments.
% % CPars=CTFitPars;
% CPars=mi.ctf;
% npanels=8;
% ndis=128;
% fig=5;
maxdist=10;

nm=size(m1,1);
pwexp=0; % pre-whitening function exponent.

nu=nm/npanels;
if ndis>nu
    ndis=nu;
end;
prewhts=Radius(nu).^pwexp;  % f^1 prewhitening
prewhts(nu/2+1,nu/2+1)=0;
cu1=single(CTF(nu,pixA,CPars(1)));
cu2=single(CTF(nu,pixA,CPars(2)));
filt=ifftshift(sign(cu1).*sign(cu2).*(abs(cu1.*cu2)).*prewhts);
% filt=ifftshift(prewhts);

ft2s=ComputeTiledFTs(m2,nu,0);
ft1s=ComputeTiledFTs(m1,nu,0);
[n1 n2 ntx nty]=size(ft2s);

ind=0;
ctu=nu/2+1;
ccs=zeros(ndis,ndis,ntx,nty);
dxs=zeros(ntx,nty);
dys=dxs;
for j=1:nty
    for i=1:ntx
        ind=ind+1;
        cc=ifftshift(real(ifftn(ft2s(:,:,i,j)...
            .*conj(ft1s(:,:,i,j)).*filt)));
        [mx x3 y3]=max2di(cc);
        dxs(i,j)=x3-ctu;
        dys(i,j)=y3-ctu;
        ccs(:,:,i,j)=imscale(Crop(cc,ndis));
    end;
end;

if fig>0
    figure(fig);
        SetGrayscale;

    distances=sqrt((dxs.^2)+(dys.^2));
    dxs(distances>maxdist)=0;
    dys(distances>maxdist)=0;
    distances=sqrt((dxs.^2)+(dys.^2));

    subplot(1,3,1);
    set(gca,'XTick',0.5:npanels+0.5);
    set(gca,'YTick',0.5:npanels+0.5);
    set(gca,'GridLineStyle','-');
    grid on;
    quiver(dxs',dys');
    
    xlabel(['Longest arrow is ' num2str(max(distances(:)*pixA)) ' A']);
    axis tight
    
subplot(1,3,2);
%     figure(fig+1);
    imacs(ImageArray(ccs));
    axis off;
%     title(name1,'interpreter','none');
    subplot(1,3,3);
    imacs(BinImage(m2c,4));
    axis off;
    
    
end
