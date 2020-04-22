% GratingDefocusMovie.m

%% Make a video of the defocus series

vlSetDarkGraphics(14);
% vlSet1080Figure;

vName='GratingsDiffr3';
makeVideo=1;

n=1024;
nfdis=[600 256]; % diffraction image size
iNdis=3*n/8; % defocused image size
nw=3*n/4; % wave display size
lw=0.5;
dirString={'under'; 'over'};
iDir=1;
nc=4; % no. of columns in plot
a1=1;
a2=1;
Cs=0; %%%%%
lambda=EWavelength(300);
showDiffractionPlane=1;


% Create the double grating
periods=[5 3];
edge=10;
vEdge=20;
nCycles=10;
boxHalf=round(0.2*n/2);
vBoxHalf=round(0.2*n/2);
vBoxShift=round(0.7*vBoxHalf);
vBoxShift=0; %%%%%

boxSize=2*boxHalf+1;
pixA=(periods(1)*nCycles)/(boxHalf*2);
xs=(-n/2:n/2-1)'*pixA;
rawGrating1=a1*cos(xs*2*pi/periods(1));
rawGrating2=a2*cos(xs*2*pi/periods(2));
msk1=fuzzymask(n,1,boxHalf/2+periods(1)/8,edge);
msk2=fuzzymask(n,1,boxHalf/2+periods(1)/8,edge);
mskv1=fuzzymask(n,1,vBoxHalf,vEdge,n/2+1-vBoxShift);
mskv2=fuzzymask(n,1,vBoxHalf,vEdge,n/2+1+vBoxShift);
grating1=rawGrating1.*msk1;
grating2=rawGrating2.*msk2;
squareGrating=grating1*mskv1'+grating2*mskv2'; % outer product

ctPars.pixA=pixA;
ctPars.lambda=lambda;
ctPars.B=0;
ctPars.alpha=.02;
ctPars.deltadef=0;
ctPars.theta=0;
ctPars.Cs=Cs;
fProj=fftn(ifftshift(squareGrating));

freqs=RadiusNorm(n)/pixA;
xFreqs=(-n/2:n/2-1)/(n*pixA);

% --------------Set up the defocus steps
dstep=.0008;
dstepX=.0002;
dmin=-.1;
dmax=.5;
iDir==1 % underfocus
ds=dmin:dstep:dmax; % Defocus values for display
dsx=dmin:dstepX:dmax;


object=squareGrating;
lineY=n/2+1;
plotLimit=1+a2;

if makeVideo
    v=VideoWriter([vName '.mp4'],'MPEG-4');
    v.FrameRate=15;
    open(v);
    disp(['Making movie file: ' vName '.mp4']);
end;

figure(1);
clf;

cxFreqs=Crop(xFreqs,nfdis(1));
cyFreqs=Crop(xFreqs,nfdis(2));
fDisc=.5;
rDisc=pixA*fDisc*n;
disc=repmat(fuzzymask(nfdis,2,rDisc,rDisc*.01)',1,1,3);
discBrightness=0.3;

if showDiffractionPlane
    subplot(121)
    %     Diffraction plane image
    uWaveSigma=3;
    uWaveAmp=n^2;
    uWave=uWaveAmp/(2*pi*uWaveSigma^2)*Gaussian(n,2,uWaveSigma);
    ctPars.defocus=-ds(1);
    cx=ContrastTransfer(freqs, ctPars, 1);
    cxImage=fftshift(fProj).*cx+uWave;
    
    fq=fctr(n);
    cxPars.x=cxFreqs;
    cxPars.y=cxFreqs;
    cxPars.scl=1e-4;
    [fdScl,rgb]=imacx2(Crop(cxImage,nfdis),1,cxPars);
    rgb1=rgb.*disc+discBrightness*(1-disc);
    image(cxFreqs,cyFreqs,rgb1);
     set(gca,'XTick',[]);
    xlabel('Spatial frequency, Å^{-1}');
end;
    % draw the static left side of the figure
    mysubplot(2,nc,1); % draw the original object
    objectC=Crop(object,iNdis);
    [objScaled,mulr,addr]=imscale(objectC,256,0);
    imaga(objectC*mulr+addr);
    axis off
if ~showDiffractionPlane
    hold on;
    plot([1 n],[lineY lineY],'-','color',[0.2 0.2 0.9],'linewidth',lw);
    hold off;

    mysubplot(4,nc,2*nc+1);
    plot(object(:,lineY));
    axis([0 inf -plotLimit plotLimit]);
end;

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
chi=-dsx'*1e4*pi*ctPars.lambda./(periods.^2)+pi/2*Cs*1e7*ctPars.lambda^4;
cFunction=sin(chi);
mysubplot(1,nc,nc);
plot(cFunction(:,1),dsx,'k-','linewidth',1.5);
hold on
plot(cFunction(:,2),dsx,'b-','linewidth',1);
hold off
axis ij

% Compute the waves ----------------------------
doComplex=1;
grating1d=sect(squareGrating);
fgrating1d=fft(ifftshift(grating1d));
wpars=struct;
wpars.x=(-nw/2:nw/2-1)*pixA;
wpars.y=ds(end:-1:1);
wpars.y=ds;
wpars.scl=0.8;
nds=numel(ds);
waves=zeros(nw,nds);
for iy=1:nds
    d=ds(iy);
    ctPars.defocus=d;
    c1=ContrastTransfer(xFreqs,ctPars);
    imgy=fftshift((ifft(fgrating1d'.*fftshift(c1))));
    waves(:,iy)=Crop(imgy,nw);
end;

% Draw the first waves
mysubplot(1,nc,nc);
cxScale=imacx2(waves,1,wpars);
% axis ij


uWaveSigma=4;
uWaveAmp=n^2*2;
uWave=uWaveAmp/(2*pi*uWaveSigma^2)*Gaussian(n,2,uWaveSigma);


for d=ds
    if showDiffractionPlane
        mysubplot(223)
        %     Diffraction plane image
%         uWaveSigma=4;
%         uWaveAmp=n^2*2;
%         uWave=uWaveAmp/(2*pi*uWaveSigma^2)*Gaussian(n,2,uWaveSigma);
        ctPars.defocus=-d;
        cx=ContrastTransfer(freqs, ctPars, 1);
        cxImage=fftshift(fProj).*cx+uWave;
        
%         fq=fctr(n);
%         cxDisImg=Crop(cxImage,n/2);
        cxPars.x=cxFreqs;
        cxPars.y=cxFreqs;
        [~,rgb]=imacx2(Crop(cxImage,nfdis),1,cxPars);
        rgb1=rgb.*disc+0.2*(1-disc);
        image(cxFreqs,cyFreqs,rgb1);
        set(gca,'XTick',[]);
    end;
        defString=dirString{iDir};
        if abs(d+.1)<.001
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
        %     title();
        text(10,20,str,'fontsize',18);
        axis off
        
      if ~showDiffractionPlane 
        %     central slice of defocused image
        mysubplot(4,nc,2*nc+2)
        plot(img(:,lineY));
        axis([0 inf -plotLimit plotLimit]);
    end;
    
    %     Sine function with marker showing our defocus level
    mysubplot(1,nc,nc-1);
    plot(cFunction(:,1),dsx,'w-','linewidth',1.5);
    hold on
    plot(cFunction(round((d-dmin)/dstepX+1),1),d,'wo','markersize',mks);
    %     plot(cFunction(:,2),dsx,'b-','linewidth',1);
    %     plot(cFunction(round((d-dmin)/dstepX+1),2),d,'bo','markersize',mks);
    hold off
    set(gca,'YTickLabel',[]);
    %     axis ij

    
    %     show complex diffracted waves
    mysubplot(1,nc,nc);
    
    cxScale=imacx2(waves,1,wpars);
    hold on;
    plot([-200 200],-ctPars.defocus*[1 1], ...
        '-','color',[1 1 1],'linewidth',lw);
    hold off;
    
    
    
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

close(v);
