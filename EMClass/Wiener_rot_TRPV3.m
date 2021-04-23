% RotSearchTRPV1.m
% Show the logP and P(\psi) for inplane rotations. Also compute it at low
% resolution to demonstrate Branch and Bound strategy. For CryoEM class on
% ML reconstruction.
% fs 22.Apr.21

doWriteFigs=0;
doLowRes=1;
fLow=.05;
rotAngle=45;
rotAngle=0
cropFactor=2;
psiAngle=180;
rng('default');


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
n=size(m0,1);
ndis=n/2;
figure(4);

kV=300;
def=1;
Cs=2;
B=100/4;
% B=0;
alpha=.05;
deltadef=0;
theta=pi/4;

m1=Crop(m0,ndis);
mc=m1/std(m1(:));
    c=ifftshift(-CTF(ndis,s.pixA,EWavelength(kV),def,Cs,B,alpha,deltadef,theta));

% make refs
rAngs=0:2:360;
nRefs=numel(rAngs);
ref0=real(ifftn(fftn(m1).*c));
refs=rsRotateImage(ref0,rAngs);
refsLow=SharpFilt(refs,fLow,fLow/10,1);
%
for irep=1:1

subplot(1,2,1);
refStep=3;
ndisC=.75*ndis;
if ~doLowRes
imagsarray(Crop(refs(:,:,1:refStep:end),ndisC,1),12);
axis off;
xlabel('References');
end;
%
sigma=25;
rImgAng=120;
    rImg0=rsRotateImage(mc,rImgAng);
    N1=(randn(ndis)*sigma);
    rImg0=real(ifftn(fftn(rImg0).*c));
    rImg=rImg0+N1;
    fc=.05;
if doLowRes
    subplot(4,4,3);
    imags(SharpFilt(rImg,fLow));
    axis off;
    subplot(4,4,4);
    imags(SharpFilt(rImg0,fLow));
    axis off;
else
    subplot(4,4,2);
    imags(GaussFilt(rImg,fc));
    axis off;
    subplot(4,4,4);
    imags(GaussFilt(rImg0,fc));
    axis off;
end;
drawnow;

%
nP=ndis^2;
rVec=reshape(refs,ndis^2,nRefs);
CC=(rImg(:)'*rVec)';
rVecLow=reshape(refsLow,ndis^2,nRefs);
CCLow=(rImg(:)'*rVecLow)';
rSq=rImg(:)'*rImg(:);

SS=rSq-2*CC+diag(rVec'*rVec);
minSS=min(SS);
% logP=-nP*log(sqrt(2*pi)*sigma)-(SS-minSS)/(2*sigma);
logP=-(SS-minSS)/(2*sigma);

maxCCLow=max(CCLow);
if ~doLowRes
    maxCCLow=NaN;
end;
logPLow=(CCLow-maxCCLow)/(4*sigma);

subplot(3,2,4);
plot([logP logPLow])
ylabel('log P(\alpha)');
subplot(3,2,6);
plot([exp(logP) exp(logPLow)]);
hold on;
plot(180-rImgAng,0,'y+');
hold off;
ylabel('P(\alpha)');
xlabel('Angle');
%%
if doLowRes
    bound=400;
    bottom=-5000;
    bottom2=-4500;
    
    subplot(3,2,1);
    plot(rAngs,[logPLow logPLow+bound])
    hold on;
    plot(rAngs,0*logPLow,'y--','linewidth',0.5);
%     plot(180-rImgAng,0,'y+');
    hold off;
    ylabel('Log(P_{low})');
    xlabel('Angle \psi');
%         legend('Log(P_{low})','Bound on Log(P)',' ','textcolor','w','fontsize',12);

    axis([0 360 bottom bound]);
    
    subplot(323);
    bottom2=-1000;
    cla;
    bins=0*rAngs+bottom2;
    bins(logPLow+bound<0)=NaN;
    bar(rAngs,bins,1,'facecolor',[.45 .22 .2]);
    hold on;
    plot(rAngs,logP,'w-');
%     zx=(rAngs);
%     zz=zx*0;
%     zz(2:2:end)=bottom;
%     zz(logPLow+bound<0)=NaN;
%     plot(zx,zz,'-','color',[.6 .4 .3]);
%     hold on;
%     plot(rAngs,logP);
     hold off;
     axis([0 360 bottom2 0]);
         ylabel('Log(P_{full})');
    xlabel('Angle \psi');

%      legend('Domain','LogP_{full}','textcolor','w');
%    
     subplot(325);
%     zz=zeros(nRefs,1);
%     zz(2:2:end)=bottom;
%     zz(logPLow<-bound)=NaN;
     plot(rAngs,exp(logP),'w-');
     axis([0 360 0 1]);
end;
ylabel('P(\alpha)');
    xlabel('Angle \psi');

% axis([0 360 -1000 0.3]);
end;


drawnow;




fname=sprintf('RotTrp%02d.jpg',irep);
fname
if doWriteFigs
print(fname,'-djpeg');
end;

return






%% -----Code left over from previous TRPV demos.










cxExponent=.4;
phaseFlip=0;
kWiener=0;
sigma=40;
showWienerCTF=kWiener>0;
color1D=[.5 .6 1];

m1=m0+sigma*randn(size(m0));
n=size(m1,1);
dFreq=1/(n*s.pixA);
kV=300;
def=1;
Cs=2;
B=100/4;
% B=0;
alpha=.05;
deltadef=0;
theta=pi/4;
c=CTF(n,s.pixA,EWavelength(kV),def,Cs,B,alpha,deltadef,theta);
if phaseFlip
    c=abs(c);
end;



% ndis=n/2;
% ndis=7*n/8;
ndis=n/2;
ndCtf=n;

nfdis=7*n/8;
mul=1/1.77;
% mul=.2;
add=130;
% sigma=0;

% set(gcf,'color',.3*[1 1 1]);
% vlSet1080Figure
vlSetDarkGraphics
figure(1);
% vlSet1080Figure;

% Object
mysubplot(2,4,1);
xs=(-ndis/2:ndis/2-1)'*s.pixA;
imags(xs,xs,Crop(m1,ndis));
text(xs(1),xs(1),'Projection',...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
% set(gca,'XTickLabel',[]);
set(gca,'XTick',[]); % Get rid of X ticks and labels
ylabel('y, Å');

mysubplot(2,4,5);
ax4=gca;

% ----------FT of object--------------
nfdis=ndis;
fRadius=s.pixA*n*0.27;
fmsk=fuzzymask(n,2,fRadius,fRadius*.2);
fm1=fmsk.*fftshift(fftn(ifftshift(m1)));
freqs=Crop((-n/2:n/2-1)'/(s.pixA*n),nfdis);
% freqs=(-ndis/2:ndis/2-1)'/(s.pixA*ndis);
pars.x=freqs;
pars.y=freqs;
imacx2(Crop(fm1,nfdis),cxExponent,pars);
text(min(freqs),min(freqs),' FT of Projection',...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
% set(gca,'XTickLabel',[]);
set(gca,'XTick',[]);
ylabel('v, Å^{-1}');

% ------Wiener filter-------

if kWiener
    cF=c./(kWiener+c.^2);
else
    cF=c;
end;


% ---------CTF---------
mysubplot(2,4,6);
if showWienerCTF
    modC=Crop(cF.*c,ndis);
    imacx2(Crop(cF.*c,ndis));
else
    modC=Crop(c,ndis);
    imacx2(Crop(c,ndis));
end;
imacx2(modC);
n1Ctf=ndis/2;
hold on;
plot(n1Ctf+[1 n1Ctf],ndis/2+1*[1 1],'w-','linewidth',1);
hold off;

text(0,0,['CTF, ' num2str(def) '\mum'],...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;

% 1D plot of CTF
mysubplot(2,4,7);
ctf1=sectr(modC);
freqs=dFreq*(1:n1Ctf);
plot(freqs,ctf1,'color',color1D);
hold on;
plot(freqs,0*ctf1,'w-','linewidth',1);
hold off;
axis([0 freqs(end) -1.05 1.05]);
text(0,-1.05,['CTF'],...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;

% ------------PSF--------------
mysubplot(2,4,2);
if showWienerCTF
    psf=Crop(fftshift(real(ifftn(ifftshift(cF.*c)))),ndis);
else
    psf=Crop(fftshift(real(ifftn(ifftshift(c)))),ndis);
end;
%     psf(fctr(ndis),fctr(ndis))=0;
if phaseFlip
    mulR=2e4;
    addR=50;
else
    [~,mulR,addR]=imscale(psf);
end;
imaga(psf*mulR+addR);
x01d=ndis/4+1;
x11d=3*ndis/4;
%     Plot the section line
hold on;
plot([x01d x11d],ndis/2+1*[1 1],'w-','linewidth',1);
hold off;
% imags(psf);
text(0,0,'PSF',...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;

% 1D plot of psf
mysubplot(2,4,3);
xs=(-ndis/2:ndis/2-1)*s.pixA;
ys=sect(psf);
ysm=max(ys);
ys=min(ysm/2,ys);

plot(xs,ys,'color',color1D);
hold on;
plot(xs,0*ys,'w','linewidth',.5);
hold off;
ySpan=max(ys)-min(ys);
y0=min(ys)-0.1*ySpan;
y1=max(ys)+0.05*ySpan;
axis([xs(1) xs(end) y0 y1]);
text(xs(1),double(y0),['PSF'],...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;

% --------------FT of image-------------
mysubplot(2,4,8);
imacx2(Crop(fm1.*cF,ndis),cxExponent);
if kWiener
    str=['FT of image, k_w = ' num2str(kWiener)];
else
    str='FT of image';
end;
text(0,0,str,...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;


% -----------Image----------------
mysubplot(2,4,4);
mc=real(fftshift(ifftn(ifftshift(fm1.*cF))));
if kWiener
    imags(Crop(mc,ndis));
else
    %     imaga(Crop(mc,ndis)/mul+128);
    imags(Crop(mc,ndis));
end;
text(0,0,'Image',...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;
% imags(mc);
return





%% % Make the demo of phase-flipping and Wiener
figure(2);
ndis=n/2;
fs=18;
def0=1;
defSpan=1;
mul=.2;

% Estimate psd of the image
sp1=RadialPowerSpectrum(m1);
spr=sp1(max(1,min(n/2-1,round(Radius(n)))));

sigma=200;
% sigma=200;
T=4;

kWiener=sigma^2./(T*spr);
% kWiener=1; %%%

defs=[def0 def0+defSpan];
B=100;
B=0;
nr=numel(defs);
nc=5;

% For combining multiple images

figure(2);
% clf;
% set(gcf,'color',.3*[1 1 1]);
mysubplot(nr,nc,1);
szs=(-ndis/2:ndis/2-1)*s.pixA;
imags(szs,szs,Crop(m0,ndis));
xlabel('angstroms');
set(gca,'YTick',[]);
% mysubplot(3,nc,nc+1);
% imags(szs,szs,-m1);
% xlabel('angstroms');
xl=ndis*.02;
yl=ndis*.01;
xu=ndis/2;
yu=ndis*.99;
for j=1:nr
    ind=j;
    N=randn(n)*sigma;
    N1=fftn(N);
    def=defs(j);
    c=CTF(n,s.pixA,EWavelength(kV),def,Cs,B,alpha,deltadef,theta);
    mc=-fftshift(real(ifftn(ifftshift(fm1.*c))));
    mysubplot(nr,nc,(j-1)*nc+2);
    imaga(Crop(GaussFilt(-mc-N,.2)*mul*.9+add,ndis));
    axis off;
    text(xl,yl,[num2str(def) '\mum defocus'],'color','w','fontsize',fs,...
        'horizontalalignment','left','verticalalignment','bottom');
    if j==1
        text(xu,yu,'Image','color','w','fontsize',fs,...
            'horizontalalignment','center','verticalalignment','top');
    end;
    
    mcf=fftshift(real(ifftn(ifftshift(fm1.*abs(c)))));
    mysubplot(nr,nc,(j-1)*nc+3);
    imaga(Crop(GaussFilt(mcf+N,.2)*mul*.9+add,ndis));
    axis off;
    if j==1
        text(xu,yu,'Phase flipped','color','w','fontsize',fs,...
            'horizontalalignment','center','verticalalignment','top');
    end;
    
    
    %     k=(sigma^2*2e-5)+1e-5;
    k=kWiener;
    wf=(fm1.*abs(c)+N1).*abs(c)./(k+abs(c).^2);
    mwf=fftshift(real(ifftn(ifftshift(wf))));
    mysubplot(nr,nc,(j-1)*nc+4);
    %     imaga(mwf*mul/1.5+add+100/(sigma+2));
    imags(Crop(mwf,ndis));
    axis off;
    if j==1
        text(xu,yu,'Wiener filtered','color','w','fontsize',fs,...
            'horizontalalignment','center','verticalalignment','top');
        cs=c;
        mcs=mc;
    else
        cs(:,:,ind)=c; % accumulate the ctfs
        mcs(:,:,ind)=mc;
    end;
    
end;
drawnow;


% Compute a combined Wiener-filtered image
ni=100;
defs=def0+defSpan*rand(1,ni);
fimgs=zeros(n,n,ni,'single');
rimgs=zeros(n,n,ni,'single');
c2s=zeros(n,n,ni,'single');
fm1=fftshift(fftn(ifftshift(m1)));
for i=1:ni
    c=-CTF(n,s.pixA,EWavelength(kV),defs(i),Cs,B,alpha,deltadef,theta);
    N1=fftn(randn(n)*sigma);
    fimgs(:,:,i)=(fm1.*c+N1).*c;
    rimgs(:,:,i)=ifftshift(real(ifftn(fftshift(fm1.*c+N1))));
    c2s(:,:,i)=c.^2;
end;
%
fWImg=sum(fimgs,3)./(kWiener+sum(c2s,3));
estImg=fftshift(ifftn(ifftshift(fWImg)));
% estImg=fftshift(ifftn(ifftshift(fm1.*c)));
mysubplot(nr,nc,nc);
imags(Crop(estImg,ndis));
% set(gca,'YTick',[]);
% xlabel('angstroms');
axis off;
text(xu,yu,['Wiener from ' num2str(ni) ' images'],'color','w','fontsize',fs,...
    'horizontalalignment','center','verticalalignment','top');
mysubplot(nr,nc,nc*nr);
c1d=sectr(sum(c2s,3));
plot(freqs,c1d(1:n1Ctf));
set(gca,'YTick',[]);

%%
% Make a big array
doRotate=1;
naDis=40;
niDis=ndis;
mc=Crop(m1,ndis);
if doRotate
    angs=360*rand(naDis,1);
    rImgs=rsRotateImage(mc,angs(1:naDis));
else
    rImgs=repmat(mc,1,1,naDis);
end;
imagsarray(rImgs);
rnImgs=zeros(ndis,ndis,naDis,'single');
for i=1:naDis
    c=ifftshift(-CTF(ndis,s.pixA,EWavelength(kV),defs(i),Cs,B,alpha,deltadef,theta));
    N1=(randn(ndis)*sigma);
    rnImgs(:,:,i)=real(ifftn(fftn(rImgs(:,:,i)).*c))+N1;
    drawnow;
end;

figure(3);
imagsarray(rnImgs);
axis off;
    



