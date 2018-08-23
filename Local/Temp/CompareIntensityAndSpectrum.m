% CompareIntensityAndSpectrum.m

% CompareIntensityAndSignal
ds=32;      % Downsampling for intensity map
dexp=.5;

imageOffset=2;  % skip this many images

cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1');
miDir=dir('Info/');
maxIndex=50;
nr=ceil(sqrt(maxIndex*.6));
nc=ceil(maxIndex/nr);

mi=ReadMiFile(['Info/' miDir(12+imageOffset).name]);
pixA=mi.pixA;
sz=mi.imageSize;
normRadius=0.1*(pixA*sz(1))/ds;  % 10 Å radius
mskd=fuzzymask(sz/ds,2,2*normRadius,.2*normRadius);

mis=cell(maxIndex,1);
mds=zeros([sz/ds maxIndex],'single');
sps=mds;
figure(1);
clf;
colormap(jet(256));
figure(2);
clf;
SetGrayscale;

for fi=1:maxIndex-imageOffset
    ind=fi+12+imageOffset;  % index in the directory

    mi=ReadMiFile(['Info/' miDir(ind).name]);
    m0=ReadMRC([mi.imagePath mi.imageFilenames{1}]);
    md=BinImage(m0,ds);
    szd=sz/ds;
    me=mean(md(:));

    mis{fi}=mi;
    mds(szd(1),szd(2),fi)=0;  % force the array to expand
    mds(:,:,fi)=md;
    
    sp=BinImage(Crop(fftshift(abs(fftn(m0-me)).^2),sz/2),ds/2);
    spmed=median(sp(:));
    spm=sp.*(1-mskd)+spmed*mskd;
    [spNorm,mulr,addr]=imscale(spm.^dexp,120);
    spNorm=sp.^dexp*mulr+addr;
    sps(szd(1),szd(2),fi)=0;
    sps(:,:,fi)=spNorm;
    
    amds=ImageArray(mds);
    asps=ImageArray(sps);
    figure(1);
    imacs(amds);
    figure(2);
    imac(asps);
    drawnow;
end;
