% CTFDemoTRPV1.m
rotAngle=45;
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
m1=rsRotateImage(squeeze(sum(mr,1)),psiAngle); % projection

%%
cxExponent=.4;
phaseFlip=1;
kWiener=.3;
showWienerCTF=0;

n=size(m1,1);
kV=300;
def=2;
Cs=2;
B=100/4;
% B=0;
alpha=.05;
deltadef=0;
theta=pi/4;
c=-CTF(n,s.pixA,EWavelength(kV),def,Cs,B,alpha,deltadef,theta);
if phaseFlip
    c=abs(c);
end;

ndis=n/2;
ndis=7*n/8;
ndis=n;
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
mysubplot(2,3,1);
xs=(-ndis/2:ndis/2-1)'*s.pixA;
imags(xs,xs,Crop(m1,ndis));
text(xs(1),xs(1),'Projection',...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
% set(gca,'XTickLabel',[]);
set(gca,'XTick',[]); % Get rid of X ticks and labels
ylabel('y, Å');

mysubplot(2,3,4);
ax4=gca;
% FT of object
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


mysubplot(2,3,5);
% CTF
if showWienerCTF
    imacx2(Crop(cF.*c,ndis));
else
    imacx2(Crop(c,ndis));
end;
text(0,0,['CTF, ' num2str(def) '\mum'],...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;

% PSF
mysubplot(2,3,2);
imags(Crop(fftshift(real(ifftn(ifftshift(c)))),ndis));
text(0,0,'Point-spread function',...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;


% FT of image
mysubplot(2,3,6);
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


% Image
mysubplot(2,3,3);
mc=real(fftshift(ifftn(ifftshift(fm1.*cF))));
if kWiener
    imags(Crop(mc,ndis));
else
    imaga(Crop(mc,ndis)/mul+128);
end;
text(0,0,'Image',...
    'color','w','fontsize',24,...
    'horizontalalignment','left','verticalalignment','bottom');
axis off;
% imags(mc);
return
%% % Make the demo of phase-flipping and Wiener

ndis=n*3/4;
fs=18;

% Estimate psd of the image
sp1=RadialPowerSpectrum(m1);
spr=sp1(max(1,min(n/2-1,round(Radius(n)))));

sigma=20;
T=4;

kWiener=sigma^2./(T*spr);
% kWiener=1; %%%

defs=[1 3];
B=100;
nr=numel(defs);
nc=5;

figure(2);
% clf;
% set(gcf,'color',.3*[1 1 1]);
mysubplot(nr,nc,1);
szs=(-ndis/2:ndis/2-1)*s.pixA;
imags(szs,szs,Crop(m1,ndis));
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
    imaga(Crop((mc+N)*mul*.9+add,ndis));
    axis off;
    text(xl,yl,[num2str(def) '\mum defocus'],'color','w','fontsize',fs,...
        'horizontalalignment','left','verticalalignment','bottom');
    if j==1
        text(xu,yu,'Image','color','w','fontsize',fs,...
            'horizontalalignment','center','verticalalignment','top');
    end;
    
    mcf=fftshift(real(ifftn(ifftshift(fm1.*abs(c)))));
    mysubplot(nr,nc,(j-1)*nc+3);
    imaga(Crop((mcf+N)*mul*.9+add,ndis));
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

%% Compute a combined Wiener-filtered image
ni=100;
defs=1+2*rand(1,ni);
fimgs=zeros(n,n,ni,'single');
c2s=zeros(n,n,ni,'single');
fm1=fftshift(fftn(ifftshift(m1)));
for i=1:ni
    c=-CTF(n,s.pixA,EWavelength(kV),defs(i),Cs,B,alpha,deltadef,theta);
    N1=fftn(randn(n)*sigma);
    fimgs(:,:,i)=(fm1.*c+N1).*c;
    c2s(:,:,i)=c.^2;
end;

fWImg=sum(fimgs,3)./(kWiener+sum(c2s,3));
estImg=fftshift(ifftn(ifftshift(fWImg)));
% estImg=fftshift(ifftn(ifftshift(fm1.*c)));
mysubplot(nr,nc,nc*nr);
imags(Crop(estImg,ndis));
% set(gca,'YTick',[]);
% xlabel('angstroms');
axis off;
        text(xu,yu,['Wiener from ' num2str(ni) ' images'],'color','w','fontsize',fs,...
            'horizontalalignment','center','verticalalignment','top');
 
