% GratingDefocusMovie2.m

%% Make a video of the defocus series
% of a 5 Å grating
% No propagating waves, just the image.

figure(1);
clf;
vlSetDarkGraphics(0,0,1); % force a reset.
vlSet1080Figure(1,1,[1024 1024]);
bkdColor=[1 1 1];
plotColor=[1 1 1];
vlSetDarkGraphics(14,bkdColor);
% vlSet1080Figure;
% plotColor=[0 .14 .2];

cd '/Users/fred/Documents/Documents - Katz/teaching/CMP710bCryoEM/710b-Videos/Video_lectures/2.2Defocus'
vName='GratingsDiffrSideBySide';
makeVideo=1;
frameRate=15;
showWaves=0;
labelFontSize=22;
n=1024;
n=2048;
% nfdis=[600 256]; % diffraction image size
iNdis=3*n/8; % defocused image size
iNdis=n/2;
plotXs=-iNdis/2:iNdis/2-1;
nw=3*n/4; % wave display size
lw=0.5;
dirString={'under'; 'over'};
iDir=1;
nc=4; % no. of columns in plot
nc=2; % no. of columns in plot
gratingAmp=[1 1];
Cs=0; %%%%%
lambda=EWavelength(300);
showDiffractionPlane=0;


% Create a single grating
periods=5;
% % Create the double grating
% periods=[5 3];
edge=20; % tapers, in pixels
vEdge=20;
nCycles=15;
frBoxHalf=0.15;
boxHalf=round(frBoxHalf*n/2);
vBoxHalf=round(frBoxHalf*0.4*n);
% vBoxShift=round(0.7*vBoxHalf);
% vBoxShift=0; %%%%%

boxSize=2*boxHalf+1;
pixA=(periods(1)*nCycles)/(boxHalf*2);
xs=(-n/2:n/2-1)'*pixA;
ngs=numel(periods);
dcOffset=.5;
squareGrating=dcOffset*fuzzymask(n,2,boxHalf*[.8 1.2],edge);
rawGrating=zeros(n,1);
msk=zeros(n,1);
mskv=zeros(n,1);
grating=zeros(n,1);
for i=1:ngs
    
    rawGrating(:,i)=gratingAmp(i)*cos(xs*2*pi/periods(i));
    msk(:,i)=fuzzymask(n,1,boxHalf/2+periods(1)/8,edge);
    mskv(:,i)=fuzzymask(n,1,vBoxHalf,vEdge,n/2+1+vBoxShift*(2*ngs-3)); % shift is -1, 1, ...
    grating(:,i)=rawGrating(:,i).*msk(:,i);
    squareGrating=squareGrating+grating(:,i)*mskv(:,i)'; % outer product
end;

ctPars.pixA=pixA;
ctPars.lambda=lambda;
ctPars.B=0;
ctPars.alpha=.02;
ctPars.alpha=0;
ctPars.deltadef=0;
ctPars.theta=0;
ctPars.Cs=Cs;%%%
ctPars.Cs=0;
fProj=fftn(ifftshift(squareGrating));

freqs=RadiusNorm(n)/pixA;
xFreqs=(-n/2:n/2-1)/(n*pixA);

% --------------Set up the defocus steps
% dstep=.0008;
% dstepX=.0002;
% dmin=-.1;
% dmax=.5;

dstep=.003;
dmin=0;
dmax1=2;
% dstep2=.1
% dmax2=2.5;
iDir==1; % underfocus
expo=2;
ds=(dmin:dstep:(dmax1^(1/expo))).^expo; % Defocus values for display
ds=[ds dmax1];
% dsx=dmin:dstepX:dmax;


object=squareGrating;
lineY=n/2+1;
lineY=iNdis/2+1;
plotLimit=1+gratingAmp(2);
plotLimit=.7+gratingAmp(2);
plotYLimits=[-plotLimit plotLimit];

if makeVideo
    v=VideoWriter([vName '.mp4'],'MPEG-4');
    v.FrameRate=frameRate;
    open(v);
    disp(['Making movie file: ' vName '.mp4']);

end;

cxFreqs=Crop(xFreqs,nfdis(1));
cyFreqs=Crop(xFreqs,nfdis(2));
fDisc=.5;
rDisc=pixA*fDisc*n;
disc=repmat(fuzzymask(nfdis,2,rDisc,rDisc*.01)',1,1,3);
discBrightness=0.3;

% if showDiffractionPlane
%     subplot(121)
%     %     Diffraction plane image
%     uWaveSigma=3;
%     uWaveAmp=n^2;
%     uWave=uWaveAmp/(2*pi*uWaveSigma^2)*Gaussian(n,2,uWaveSigma);
%     ctPars.defocus=-ds(1);
%     cx=ContrastTransfer(freqs, ctPars, 1);
%     cxImage=fftshift(fProj).*cx+uWave;
%     
%     fq=fctr(n);
%     cxPars.x=cxFreqs;
%     cxPars.y=cxFreqs;
%     cxPars.scl=1e-4;
%     [fdScl,rgb]=imacx2(Crop(cxImage,nfdis),1,cxPars);
%     rgb1=rgb.*disc+discBrightness*(1-disc);
%     image(cxFreqs,cyFreqs,rgb1);
%     set(gca,'XTick',[]);
%     xlabel('Spatial frequency, Å^{-1}');
% end;
% draw the static left side of the figure
mysubplot(2,nc,1); % draw the original object
objectC=Crop(object,iNdis);
[objScaled,mulr,addr]=imscale(objectC,256,0);
imaga(objectC*mulr+addr);
axis off
% if ~showDiffractionPlane
    hold on;
    plot([1 n],[lineY lineY],'-','color',[0.2 0.2 0.9],'linewidth',lw);
    hold off;
    
    mysubplot(4,nc,2*nc+1);
    plot(plotXs*pixA,objectC(:,lineY));
    set(gca,'Color',plotColor);
    axis([-inf inf plotYLimits]);
        xlabel('X coordinate, Å','Color','k','FontSize',labelFontSize);
% end;

% if makeVideo
%     f=getframe(gcf);
%     writeVideo(v,f);
%     imwrite(f.cdata,[vName '1.jpg']);
% end;

addr=180;
% CTF(n,Pars)
% a struct Pars containing each of these fields is passed:
% pixA (optional); lambda, defocus, Cs, B, alpha, deltadef, theta.

mks=15;

% chi=-1e4*lambda.*defocus.*s2+Cs.*lambda.^3.*5e6.*s2.^2-alpha/pi
% chi=-dsx'*1e4*pi*ctPars.lambda./(periods.^2)+pi/2*Cs*1e7*ctPars.lambda^4;
% cFunction=sin(chi);
% if showWaves
% mysubplot(1,nc,nc);
% plot(cFunction(:,1),dsx,'k-','linewidth',1.5);
% hold on
% for i=2:ngs
%     plot(cFunction(:,2),dsx,'b-','linewidth',1);
% hold off
% end;
% axis ij
% end;

% % Compute the waves ----------------------------
% doComplex=1;
% grating1d=sect(squareGrating);
% fgrating1d=fft(ifftshift(grating1d));
% wpars=struct;
% wpars.x=(-nw/2:nw/2-1)*pixA;
% wpars.y=ds(end:-1:1);
% wpars.y=ds;
% wpars.scl=0.8;
% nds=numel(ds);
% waves=zeros(nw,nds);
% for iy=1:nds
%     d=ds(iy);
%     ctPars.defocus=d;
%     c1=ContrastTransfer(xFreqs,ctPars);
%     imgy=fftshift((ifft(fgrating1d'.*fftshift(c1))));
%     waves(:,iy)=Crop(imgy,nw);
% end;
% 
% % Draw the first waves
% if showWaves
%     mysubplot(1,nc,nc);
%     cxScale=imacx2(waves,1,wpars);
%     % axis ij
% end;

% uWaveSigma=4;
% uWaveAmp=n^2*2;
% uWave=uWaveAmp/(2*pi*uWaveSigma^2)*Gaussian(n,2,uWaveSigma);


for d=ds
%     if showDiffractionPlane
%         mysubplot(223)
%         %     Diffraction plane image
%         %         uWaveSigma=4;
%         %         uWaveAmp=n^2*2;
%         %         uWave=uWaveAmp/(2*pi*uWaveSigma^2)*Gaussian(n,2,uWaveSigma);
%         ctPars.defocus=-d;
%         cx=ContrastTransfer(freqs, ctPars, 1);
%         cxImage=fftshift(fProj).*cx+uWave;
%         
%         %         fq=fctr(n);
%         %         cxDisImg=Crop(cxImage,n/2);
%         cxPars.x=cxFreqs;
%         cxPars.y=cxFreqs;
%         [~,rgb]=imacx2(Crop(cxImage,nfdis),1,cxPars);
%         rgb1=rgb.*disc+0.2*(1-disc);
%         image(cxFreqs,cyFreqs,rgb1);
%         set(gca,'XTick',[]);
%     end;
    defString=dirString{iDir};
    if abs(d+.1)<.001 && markScherzer
        ctPars.defocus=.088;
        defString='(Scherzer)';
    else
        ctPars.defocus=-d+eps;
    end;
    c=CTF(n,ctPars);
    fimg=fProj.*fftshift(c);
    img=fftshift(real(ifftn(fimg)));
    imgC=Crop(img,iNdis);
    
    
    mysubplot(2,nc,2); % show the defocused image
    addr1=addr-50;
    imaga(mulr*imgC+addr1);
    hold on;
    plot([1 n],[lineY lineY],'-','color',[0.2 0.2 0.9],'linewidth',lw);
    hold off;
    str=['focus ' sprintf('%04.3f',ctPars.defocus) '\mum ' defString];
    str=['defocus: ' sprintf('%04.2f',-ctPars.defocus) '\mum ' defString];
    %     title();
%     text(10,20,str,'fontsize',18);
    text(50,50,str,'fontsize',24);
    axis off
    
%     if ~showDiffractionPlane
        %     central slice of defocused image
        mysubplot(4,nc,2*nc+2)
        plot(plotXs*pixA,imgC(:,lineY));
        set(gca,'Color',plotColor);
        axis([-inf inf plotYLimits]);
        xlabel('X coordinate, Å','Color','k','FontSize',labelFontSize);
%     end;
    
%     if showWaves
%         %     Sine function with marker showing our defocus level
%         mysubplot(1,nc,nc-1);
%         plot(cFunction(:,1),dsx,'w-','linewidth',1.5);
%         hold on
%         plot(cFunction(round((d-dmin)/dstepX+1),1),d,'wo','markersize',mks);
%         %     plot(cFunction(:,2),dsx,'b-','linewidth',1);
%         %     plot(cFunction(round((d-dmin)/dstepX+1),2),d,'bo','markersize',mks);
%         hold off
%         set(gca,'YTickLabel',[]);
%         %     axis ij
%         
%         
%         %     show complex diffracted waves
%         mysubplot(1,nc,nc);
%         
%         cxScale=imacx2(waves,1,wpars);
%         hold on;
%         plot([-200 200],-ctPars.defocus*[1 1], ...
%             '-','color',[1 1 1],'linewidth',lw);
%         hold off;
%     end;
%     
    
    %     nc0=floor(nc/2);
    %     mysubplot(2,nc0,nc0+1);
    %     imacx2(fftshift(fimg),.2);
    drawnow;
    if makeVideo
        f=getframe(gcf);
        writeVideo(v,f);
        %     imwrite(f.cdata,[vName '1.jpg']);
    end;
    %     pause(0.1);
end;

if makeVideo
close(v);
end;