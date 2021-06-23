% WienerTRPV3.m
% Make a movie of accumulating noisy images

cd ~/aEMCodeRepository/EMClass/
%%
% vlSetDarkGraphics;
% vlSet1080Figure;

vName='WienerTRPV3r2constK';
frameRate=2;
makeVideo=0;
cropFactor=2;
load emd_5778_projs.mat
m0=Crop(mr0,512);  % side view, padded
pixA=s.pixA;

n=size(m0,1);
dFreq=1/(n*s.pixA);
kV=300;
% def=1;
Cs=2;
B=100/4;
B=50/4;
% B=0;
alpha=.05;
deltadef=0;
theta=pi/4;
c=CTF(n,pixA,EWavelength(kV),def,Cs,B,alpha,deltadef,theta);

% Compute a combined Wiener-filtered image
constWiener=100; % Don't use the SSNR, use this value.
constWiener=0; % use the SSNR
ni=3e3;
dCrop=8/3;
showCis=1;
fc=.4;
ndis=n/dCrop;
nr=2;
nc=3;
def0=1;
defSpan=.5;
xu=ndis/2;
yu=ndis*.99;
fs=18;
plotColor=[.8 .8 1];
plotColor2=[.3 .3 .6];
plotColor3=[.6 .6 .6];

origColorOrder=colororder;

reps=[10 15 20 30 50 70];
niDisVals=[1:9 reps 10*reps 100*reps 1000*reps 1e4*reps 1e6];
nNiDis=numel(niDisVals);
freqs1=(1:n/2)'/(pixA*n);

d1=fuzzymask(n,2,60);  % tight mask
q=d1.*m0;
mVar=(q(:)'*q(:))/(d1(:)'*d1(:));
m1=m0/sqrt(mVar); % normalized to unit power.

% ---show the original projection----
subplot(nr,nc,1);
xs=ndis*pixA/2*[-1 1];
imags(xs,xs,Crop(m1,ndis));
xlabel('X coordinate, Å');
ylabel('Y coordinate, Å');
drawnow;

% Estimate psd of the image
sp1=RadialPowerSpectrum(d1.*m1);
sp1(end)=sp1(end-1);
sp1=Smooth1D(sp1,20);
spr=ToRect(sp1,n);
% subplot(nr,nc,2);
% semilogy(sp1);

sigma=.6;

T=2;

SSNR1D=sp1/sigma^2;
kWiener1D=1./(T*dCrop^2*SSNR1D);
if constWiener>0
    kWiener1D=0*sp1+constWiener;
end;

% --Show SSNR--
subplot(nr,nc,nc+1);
ssnr1=1./kWiener1D;
semilogy(freqs1,ssnr1);
axis([0 inf min(ssnr1)*.7 max(ssnr1)*1.2 ]);
ylabel('Spectral SNR');
xlabel('Spatial frequency');
drawnow;

kWiener=ToRect(kWiener1D,n);

niDisInd=1;

defs=def0+defSpan*rand(1,ni);
fimgAcc=zeros(n,n,'single');
rimgAcc=zeros(n,n,'single');
c2Acc=zeros(n,n,'single');
allFRCs=0.5*ones(n/2,numel(niDisVals));
Cis=zeros(n/2,ni,'single');
fm1=fftshift(fftn(ifftshift(m1)));

if makeVideo
    v=VideoWriter([vName '.mp4'],'MPEG-4');
    v.FrameRate=frameRate;
    open(v);
    disp(['Making movie file: ' vName '.mp4']);
end;


for i=1:ni
    c=-CTF(n,s.pixA,EWavelength(kV),defs(i),Cs,B,alpha,deltadef,theta);
    N1=fftn(randn(n)*sigma);
    fimg=fm1.*c+N1;
    fimgAcc=fimgAcc+c.*fimg;
    c2Acc=c2Acc+c.^2;
    Cis(:,i)=sectr(c).^2;
    if i==1
        rimg=GaussFilt(-ifftshift(real(ifftn(fftshift(fimg)))),fc);
        [~,mulr,addr]=imscale(rimg,256);
    end;
    if i>=niDisVals(niDisInd) % do the display update.

        %         ----Show the raw image----
        rimg=GaussFilt(-ifftshift(real(ifftn(fftshift(fimg)))),fc);
        subplot(nr,nc,2);
        imaga(addr+mulr*Crop(rimg,ndis));
        axis off;
        
        %         ---Show the estimated image---        
        fWImg=fimgAcc./(kWiener+c2Acc);
        estImg=fftshift(ifftn(ifftshift(fWImg)));
        % estImg=fftshift(ifftn(ifftshift(fm1.*c)));        
        subplot(nr,nc,nc);
        imags(Crop(estImg,ndis));
        % set(gca,'YTick',[]);
        % xlabel('angstroms');
        axis off;
        title(['N= ' num2str(i) ' images'],'color','w');
        %         text(xu,yu,['Wiener from ' num2str(i) ' images'],'color','w','fontsize',fs,...
        %             'horizontalalignment','center','verticalalignment','top');
        
        % ----plot the denominator terms----
        subplot(nr,nc,nc+2);
        c1d=sectr(c2Acc);
        %         plot(freqs1,c1d(1:n1Ctf)/i);
        kw=kWiener1D/i;
        %         kw(kw>.5)=nan; % Clip the KWiener
        plot(freqs1,[c1d/i kw]);
        axis([0 inf 0 1 ]);
        xlabel('Spatial frequency');
        ylabel('Denominator terms /N');
        legend('Sum of CTF^2','k_{Wiener}','textcolor','w')
        %         set(gca,'YTick',[]);

%         ----show the sum of Ci^2----
        if showCis && i>1
            nf=n/4;
            iLim=min(300,i);
            subplot(nr,nc,nc+1);
            plot(freqs1(1:n/4),.5+.5*Cis(1:n/4,1:iLim));
            hold on;
            plot(freqs(1:nf),.5*c1d(1:nf)/i,'color',plotColor);
            hold off;
            axis([0 inf 0 1]);
            set(gca,'ytick',[]);
            ylabel('Sum of CTF^2');
            xlabel('Spatial frequency');
   text(freqs1(nf)*.95,.95,['\Deltaz=' num2str(defs(i),3)],'color','w','fontsize',fs,...
              'horizontalalignment','right','verticalalignment','top');
%                     title(['\Delta z = ' num2str(defs(i),3) ' \mum'],'color','w');

        end;
        
        %         ----plot the FRCs----
        subplot(nr,nc,nc*nr)
        frc=FRC(d1.*m0,d1.*estImg);
        smFrc=Smooth1D(frc,30);
        allFRCs(:,niDisInd)=Smooth1D(smFrc,50);
        if niDisInd>1
            plot(freqs1,allFRCs(:,1:niDisInd-1),'color',plotColor2);
        else
            plot(freqs1,nan+0*smFrc);
        end;
        hold on;
        plot(freqs1,smFrc*0+.5,'--','color',plotColor3);
        %         plot(freqs1,[smFrc*0 smFrc*0+1],'color',plotColor3);
        plot(freqs1,smFrc,'color',plotColor);
        hold off;
        axis([0 inf 0 1]);
        ylabel('Fourier ring correlation');
        xlabel('Spatial frequency');
        title(['N= ' num2str(i)],'color','w');
        
        drawnow;
        if makeVideo
            f=getframe(gcf);
            writeVideo(v,f);
        end;
        niDisInd=niDisInd+1;
        %              pause
    end;
end;

if makeVideo
    close(v);
    disp('done.');
end;

%

