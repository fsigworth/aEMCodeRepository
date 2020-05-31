% Figs for intro lectures

% TRPV1 map for displaying defocus contrast
cd('/Users/fred/Documents/Documents - Katz/papers/Microscopy review/Figs/EMPIAR 10005/EMD-5778/map');
[m,s]=ReadMRC('emd_5778.map');
map=Downsample(m,128);
rotMapd=(ERotate3(map,[0 pi pi/2]));
ShowSections(rotMapd);
% WriteMRC(rotMapd,s.pixA*2,'RotTRPV128.map');

%%
rotMap=ERotate3(m,[0 pi pi/2]);
% WriteMRC(rotMap,s.pixA,'RotTRPV256.map');
%%
% dfScherzer=-1.2*sqrt(Cs*lambda)
Cs=2.7e-3;
lambda=.02e-10;
dfScherzer=1.2*sqrt(Cs*lambda)
% =88 nm.
% resolution:
% u_{{\text{res}}}({\text{Scherzer}})=0.6\lambda ^{{3/4}}C_{s}^{{1/4}},
uRes=0.6*lambda^.75*Cs^.25
%% map, power spectrum, zero defocus
invB=20;
ndis=192;
n=256;

[x,y,z]=ndgrid(-n/2:n/2-1);
f3=sqrt(x.^2+y.^2+z.^2)/(n*s.pixA);
% filt3=exp(50*f3.^2-(f3/.3).^6).*fuzzymask(n,3,.35*n);;
filt3=exp(invB*f3.^2).*fuzzymask(n,3,.4*n,.1*n);
% filt3=1;
rotMapSharp=real(ifftn(fftn(rotMap).*fftshift(filt3)));
sp1=RadialPowerSpectrum(Crop(rotMapSharp,128));
fs=(0:63)'/(128*s.pixA);
filt=exp(invB*fs.^2-(fs/.3).^6);
figure(2);
mysubplot(2,2,2);
semilogy(fs,[sp1 sp1]);
drawnow;

%% Make the masked map for ctf demo
lpMap=1*GaussFilt(rotMapSharp,.005);
xProj=squeeze(sum(-rotMapSharp-lpMap,1));
xProj=SharpFilt(xProj,.4);
% Make a tighter mask
msk=fuzzymask(n,2,n/4,n/8);
q=(msk.*(xProj-100));
msk2=GaussFilt(GaussFilt(q<-100,.1)>.1,.02);
xProjC=Crop(xProj.*msk2,ndis);
clf;
imags(xProjC);

%% Make a video of the defocus series
showTrp=1;
lineY=ndis/2+10;
lw=0.5;
plotScale=1/800;
dirString={'over'; '' ; 'under'};
makeVideo=0;

% Create a double grating: GratingDefocusMovie does this better!
nd=1024;
per=5; % Å
per2=2;
edge=10;
vEdge=20;
a2=1;
nCycles=5;
boxHalf=round(0.1*nd/2);
vBoxHalf=round(0.15*nd/2);
vBoxShift=round(0.7*vBoxHalf);
boxSize=2*boxHalf+1;
pixA=(per*nCycles)/(boxHalf*2);
xs=(-nd/2:nd/2-1)'*pixA;
rawGrating1=cos(xs*2*pi/per);
rawGrating2=a2*cos(xs*2*pi/per2);
msk=fuzzymask(nd,1,boxHalf+per/8,edge);
mskv1=fuzzymask(nd,1,vBoxHalf,vEdge,nd/2+1-vBoxShift);
mskv2=fuzzymask(nd,1,vBoxHalf,vEdge,nd/2+1+vBoxShift);
grating1=rawGrating1.*msk;
grating2=rawGrating2.*msk;
squareGrating=grating1*mskv1'+grating2*mskv2'; % outer product


if showTrp
    object=-xProjc;
    pixA=s.pixA;
    plotScale=.0013;
    ndis=192;
else
    object=squareGrating;
    lineY=ndis/2+1;
    plotScale=1/(1+a2);
    ndis=nd;
end;

vName='Defocus';
if makeVideo
    v=VideoWriter([vName '.avi']);
    v.FrameRate=15;
    open(v);
    disp(['Making movie file: ' vName '.avi']);
end;

figure(2);
clf;

n=size(rotMap,1);

% draw the static left side of the figure
mysubplot(2,2,1);
[objScaled,mulr,addr]=imscale(object,256,0);

imaga(object*mulr+addr);
axis off equal
hold on;
plot([1 ndis],[lineY lineY],'-','color',[0.2 0.2 0.9],'linewidth',lw);
hold off;
mysubplot(4,2,5);
plot(object(:,lineY)*plotScale);
axis([0 inf -.15 .85]);

if makeVideo
    f=getframe(gcf);
    writeVideo(v,f);
    imwrite(f.cdata,[vName '1.jpg']);
end;
    
addr=180;
% CTF(n,Pars)
% a struct Pars containing each of these fields is passed:
% pixA (optional); lambda, defocus, Cs, B, alpha, deltadef, theta.
ctPars.pixA=pixA;
ctPars.lambda=EWavelength(300);
ctPars.B=0;
ctPars.alpha=.02;
ctPars.deltadef=0;
ctPars.theta=0;
ctPars.Cs=2.7;
fProj=fftn(ifftshift(object));
freqs=RadiusNorm(n)/s.pixA;

for d=0:.05:.5
    ind=2+sign(-d);
    defString=dirString{ind};
    if abs(d+.1)<.001
        ctPars.defocus=.088;
        defString='(Scherzer)';
    else
        ctPars.defocus=-d+eps;
    end;
    c=CTF(ndis,ctPars);
    img=fftshift(real(ifftn(fProj.*fftshift(c))));

%     plot(sectr(c));
    mysubplot(2,2,2);
    imaga(xs,xs,mulr*img+addr-20);
hold on;
    plot(xs,0*xs,'-','color',[.5 .5 .9],'linewidth',lw);
hold off;

str=['focus ' sprintf('%04.3f',ctPars.defocus) '\mum ' defString];
%     title();
    text(10,20,str,'fontsize',18);
    axis off equal
    
    mysubplot(4,2,6)
    plot(img(:,lineY)*plotScale);
    axis([0 inf -.5 .5]);
    
    pause(0.1);
end;




