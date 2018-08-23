% VesicleModelFigure
% 
nr=2;
nc=4;

cd('/Users/fred/EMWork/Hideki/160909p/KvLipo121_2_2picking');
mi=ReadMiFile('Info/sq02_1_0004_Sep09_19.06.23mi.txt');
m=ReadMRC([mi.procPath mi.baseFilename 'ms.mrc']);
vi=34;  % vesicle index;
ndis=220;
ds=2;
nm=mi.imageSize/ds;
%%
miRound=mi;
miRound.vesicle.r(:,2:end)=0;

miConst=miRound;
miConst.vesicle.s(:,2:end)=0;

miConstShell=miConst;
miConstShell.vesicleModel=zeros(mw,1);
miConstShell.vesicleModel(mw2+1)=1;

miConstEl=mi;
miConstEl.vesicle.s(:,2:end)=0;

miFull=mi;
miFull.vesicle.s(:,:,2:end)=0;

miCorr=mi;
miCorr.vesicle.s(:,:,1)=0;
miCorr.vesicle.s(:,2:end,:)=0;

vConstShell=ExtractImage( meMakeModelVesicles(miConstShell,nm,vi,0,0),...
    round([mi.vesicle.x(vi) mi.vesicle.y(vi)]/ds)+1,ndis);
vConst=ExtractImage( meMakeModelVesicles(miConst,nm,vi,0,0),...
    round([mi.vesicle.x(vi) mi.vesicle.y(vi)]/ds)+1,ndis);
vConstEl=ExtractImage( meMakeModelVesicles(miConstEl,nm,vi,0,0),...
    round([mi.vesicle.x(vi) mi.vesicle.y(vi)]/ds)+1,ndis);
vFull=ExtractImage( meMakeModelVesicles(miFull,nm,vi,0,0),...
    round([mi.vesicle.x(vi) mi.vesicle.y(vi)]/ds)+1,ndis);
vFullC=ExtractImage( meMakeModelVesicles(miFull,nm,vi,1,0),...
    round([mi.vesicle.x(vi) mi.vesicle.y(vi)]/ds)+1,ndis);
vCorr=ExtractImage( meMakeModelVesicles(miCorr,nm,vi,0,0),...
    round([mi.vesicle.x(vi) mi.vesicle.y(vi)]/ds)+1,ndis);


%%
figure(2);

subplot(nr,nc,nc+1);
mw=numel(mi.vesicleModel);
mw2=floor(mw/2);
plot((-mw2:mw2)*mi.pixA,mi.vesicleModel);
axis([-inf inf 0 2.5]);
xlabel('Membrane normal, Å');
ylabel('Inner potential, V');
title('Membrane scattering');

xs=(1:ndis)*mi.pixA*ds;
subplot(nr,nc,1);
imags(xs,xs,-vConstShell);
xlabel('X coordinate, Å');
ylabel('Y coordinate, Å');
title('Spherical shell');

subplot(nr,nc,2);
imags(-vConst);
axis off;
title('Membrane density');

subplot(nr,nc,3);
imags(-vConstEl);
[~,mulr,addr]=imscale(-vConstEl);
axis off;
title('+variable radius');

subplot(nr,nc,4);
imags(-vFull);
axis off;
title('+variable density');

subplot(nr,nc,2*nc-1);
imags(vFullC);
title('+CTF');
axis off;

dsx=mi.imageSize(1)/size(m,1);
mx=ExtractImage(m,round([mi.vesicle.x(vi) mi.vesicle.y(vi)]/dsx)+1,ndis*ds/dsx);
subplot(nr,nc,2*nc-2);
imags(mx);
axis off;

subplot(nr,nc,2*nc);
imags(mx-Downsample(vFullC,ndis*ds/dsx));
axis off;

