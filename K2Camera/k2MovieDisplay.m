% k2MovieDisplay
% Show images and spectra for the frames of a movie.
% Make a plot showing the time of a focal jump.

pixA=1.25;
groupSize=5;

[fname, pa]=uigetfile('*.*','Select a movie files');
if isnumeric(pa) % File selection cancelled
    return
end;
cd(pa);

disp(['Reading ' fname]);
mv=ReadMRC(fname);
disp(' done.');
fmv=single(mv);

mvMean=mean(fmv(:));
n=[size(mv,1) size(mv,2)];

n=[3840 3840];

fmv=Crop(fmv,n,1,mvMean);  % pad with the mean value


nf=size(mv,3);  % number of frames
disp([num2str(nf) ' frames']);
nr=floor(sqrt(nf));
nc=nr;
if nr*nc<nf
    nc=nc+1;
end;
if nr*nc<nf
    nr=nr+1;
end;

%%
nFirstAvg=floor(nf/5);
nLastAvg=floor(nf/5);
nFirstAvg=1;
nLastAvg=1;

%  Display parameters
imBinFactor=8;
imCropFactor=4;
spBinFactor=8;
spCropFactor=2;
spExponent=.3;
spTruncVals=[0 .00];  % fraction of black and white values to truncate
figure(5);
SetGrayscale;
imac(imscale(BinImage(sum(fmv,3),imBinFactor),256,[0 .00011]));

figure(1);
clf;
SetGrayscale;
nbx=1:imBinFactor:n(1);
nby=1:imBinFactor:n(2);
ncs=round(n/imBinFactor/imCropFactor);
xs=Crop(nbx,ncs(1));
ys=Crop(nby,ncs(2));
sims=[];
for i=1:nf
    subplot(nr,nc,i);
    sim=Crop(BinImage(fmv(:,:,i),imBinFactor),round(n/imBinFactor/imCropFactor));
    imacs(xs,ys,sim);
    sims(:,:,i)=sim;
    if i==1
        title(fname,'interpreter','none');
    else
        title(i);
    end;
end;
%
firstAvg=mean(sims(:,:,1:nFirstAvg),3);
firstAvg=firstAvg-mean(firstAvg(:));
lastAvg=mean(sims(:,:,end-nLastAvg+1:end),3);
lastAvg=lastAvg-mean(lastAvg(:));
correls=zeros(nf,1);
for i=1:nf
    sim=sims(:,:,i);
    sim=sim-mean(sim(:));
    correls(i,1)=sim(:)'*firstAvg(:);
    correls(i,2)=sim(:)'*lastAvg(:);
end;
figure(4);
plot(correls,'.-','markersize',10);
drawnow;

%%
ncr=n(1)/spBinFactor/spCropFactor;  % reduced spectrum size
df=spBinFactor/(pixA*n(1));  % frequency step
f1=.03/df;  % 30A  % radius for this frequency
f2=.08/df;  % 12A;
band=fuzzymask(ncr,2,f2,.1*f2)-fuzzymask(ncr,2,f1,.5*f1);
spExponent=.3;
%
figure(2);
clf;
SetGrayscale;
ncr=round(n/spBinFactor/spCropFactor);
xs=(-ncr(1)/2:ncr(1)/2-1)/(ncr(1)*spCropFactor);
ys=(-ncr(2)/2:ncr(2)/2-1)/(ncr(2)*spCropFactor);
s1d=[];
fmfb0=zeros(ncr);
% correls=zeros(nf,1);
sims=zeros([ncr nf]);
for i=1:nf
    mf=fmv(:,:,i);
    mf=mf-mean(mf(:));
    mf=mf-mean(mf(:));
    fmf=Crop(fftshift(fftn(mf)),round(n/spCropFactor));
    fmfb=BinImage(fmf,spBinFactor);
    
    spReduced=abs(fmfb);
    sps(:,:,i)=spReduced;
    s1d(:,i)=Radial(spReduced);
%     figure(2);
%     subplot(nr,nc,i);
%     imac(xs,ys,sim);
%     if i==1
%         title(fname,'interpreter','none');
%     else
%         title(i);
%     end;
%     drawnow;
    q=fmfb.*band;
    correls(i)=real(q(:)'*fmfb0(:));
    fmfb0=fmfb;
figure(4);
plot(correls,'.-');
drawnow;
end;
%
figure(2);
for i=1:nf;
        sim=imscale(sps(:,:,i).^spExponent,255,spTruncVals);
    subplot(nr,nc,i);
    imac(xs,ys,sim);
end;
    % figure(3);
% semilogy(s1d);

%%
legend(num2str((1:nf)'));

