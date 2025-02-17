% WienerTRPV2.m
% demonstration of phase-flipping, Wiener on images.
% Version 2, plots 1D CTF as well.

cd ~/aEMCodeRepository/EMClass/
skipProjections=1;


if ~skipProjections
% Code to make projections
psiAngle=180;
[m,s]=ReadMRC('emd_5778.map');
rotAngle=0;
for iRot=0:1
    rotAngle
    if mod(rotAngle,90)~=0
        disp('Rotating...');
        mr=rsRotateImage(m,rotAngle); % rotate about z axis
        disp('done.');
    else
        mr=m;
    end;
    mr=rsRotateImage(squeeze(sum(mr,1)),psiAngle); % projection
    
    % Solvent flattening
    mr1=35+(mr+4*GaussFilt(mr,.0027));
    mr2=mr1.*fuzzymask(256,2,100,20);
    % plot(sect(m0))
    
    mTop=sum(m,3);
    mt1=35+(mTop+4*GaussFilt(mTop,.0027));
    mt2=mt1.*fuzzymask(256,2,100,20);
    mt0=Crop(mt2,512);  % top view
    % imags(mt0);
    
    if rotAngle==0
        mr0=mr2;
    else
        md0=mr2;
    end;
    rotAngle=45;
end;
mt0=mt2;

% write out X and diagonal side projections, and top view.
WriteMRC(mr0,s.pixA,'emd_5778_yz.mrc'); % all 256 pixels
WriteMRC(md0,s.pixA,'emd_5778_dz.mrc');
WriteMRC(mt0,s.pixA,'emd_5778_xy.mrc');
save('emd_5778_projs.mat','mr0','md0','mt0','s');
end;
%%
vlSetDarkGraphics;
% vlSet1080Figure;

    cropFactor=2;
    [mr0,s]=ReadMRC('emd_5778.map');
load emd_5778_projs.mat
    m0=Crop(mr0,512);  % side view, padded
    pixA=s.pixA;
    
% Which version to make:
phaseFlip=0;
kWiener=0;

sigma=40;
showWienerCTF=kWiener>0;
color1D=[.7 .8 1];
color1D2=[.4 .5 .7];
cxExponent=.4;
% color1D2=color1D;

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
c=CTF(n,pixA,EWavelength(kV),def,Cs,B,alpha,deltadef,theta);
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
imacx2(Crop(fm1,nfdis),cxExponent*.8,pars);
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
imacx2(modC); % Show the net ctf 
n1Ctf=ndis/2;
hold on;
plot(n1Ctf+[1 n1Ctf],ndis/2+1*[1 1],'w-','linewidth',1);
hold off;

text(0,0,['CTF, ' num2str(def) '\mum'],...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;

% -------1D plot of CTF------------
mysubplot(2,4,7);
ctf1=sectr(modC);
freqs=dFreq*(1:n1Ctf);
if phaseFlip
    c1=sectr(c);
    plot(freqs,c1(1:n1Ctf),'--','color',color1D2);
    hold on;
    plot(freqs,ctf1,'color',color1D);
else
    plot(freqs,ctf1,'color',color1D);
    hold on;
end;
plot(freqs,0*ctf1,'w-','linewidth',1);
hold off;
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
    addR=70;
else
    [~,mulR,addR]=imscale(psf);
end;
imaga(psf*mulR+addR); % Show the 2D psf

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
if phaseFlip
ysm=max(ys);
ys=min(ysm/4,ys);
end;

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
n=256;
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
    



