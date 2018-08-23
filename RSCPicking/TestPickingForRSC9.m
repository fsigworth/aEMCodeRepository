% function TestPickingForRSC9 %(m,mi,map,mapPixA)
% Generic inverse filter params
f0=.02;  % 4th order hp
% f0=.05;
ex=4;
f1=.0015; % Lorentizian
% f1=.004; % Lorentizian
a1=.4;

ds1=2;  % extra downsampling compared to the merged image.
simulateImage=0;
showTemplates=0;

% mapName='/Volumes/TetraData/Structures/AMPAR/3KG2map58.mrc';
mapName='~/aEMCodeRepository/Data/3KG2map58.mrc';

[fname pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
[rootPath infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;
cd(rootPath);

fileIndex=1;
miName=[infoPath fname{fileIndex}];
load(miName);
if ~isfield(mi,'vesicle.ok')
    amps=mi.vesicle.s;
    minS=0;
    maxS=inf;
    mi.vesicle.ok=(amps>minS) & (amps<maxS);
end;

imgExts={'mc.mrc','m.mrc'};
i=0;
ok=0;
while i<=numel(imgExts)&&~ok
    i=i+1;
    imgName=[mi.procPath mi.baseFilename imgExts{i}];
    if exist(imgName,'file')
        ok=1;
        origImg=ReadEMFile(imgName);
    end;
end;
if ~ok
    error(['Couldn''t open the image ' imgName]);
end;

n=size(origImg);

if ds1>1  % further downsampling
    n=n/ds1;
    m0=Downsample(origImg,n);
else
    m0=origImg;
end;

ds=mi.imageSize(1)/n(1);
pixA=mi.pixA*ds;

disp('Making model vesicles');
ves=meMakeModelVesicles(mi,n);

m1=m0-ves;

disp('Noise whitening');
H=meGetNoiseWhiteningFilter(mi,size(m1));
m=real(ifftn(fftn(m1).*ifftshift(H)));

figure(1); clf;
SetGrayscale;
imac(imscale(GaussFilt(m,.2),256,.001));
drawnow;

% Load the 3D map
disp('Loading the 3D map');
[origMap s]=ReadMRC(mapName);
mpixA=s.pixA;
nt1=size(origMap,1)*mpixA/pixA;  % final effective map size
nt=ceil(nt1/8)*8;
[map finalmag]=DownsampleGeneral(origMap,nt,mpixA/pixA);
% magnifications=[mpixA/pixA finalmag]

%% Create the list of angles for the templates
ds=2*ds1;  % We started with a 2x downsampled merged image.  Overall downsampling is 4.
membraneOffset=-24/ds;  % downsampled map by 2.
nAlpha=32; % about 10 degrees
nBeta=12;  % even is best.  Number on hemisphere.
nGamma=8;
symmetry=2;
gammaStep=360/(symmetry*nGamma);
nterms=26;
nterms=12;

[hemiAngles angleInds]=rsListHemisphereAngles(nAlpha, nBeta);
nHemiAngles=size(hemiAngles,1);
nHemi=2;  % both hemispheres

angleList=zeros(nHemi,nGamma,nHemiAngles,3);
for j=1:nGamma;
    gamma=(j-1)*gammaStep;
    for k=1:nHemiAngles
        angleList(1,j,k,:)=[hemiAngles(k,:) gamma];
        angleList(2,j,k,:)=[[0 180]-hemiAngles(k,:) gamma];
    end;
end;
nAngles=numel(angleList)/3;

% angle list is of size
% (nHemi x nGamma x nHemiAngles, 3)  where nHemi=2 is the number of hemispheres.
% Note that the beta angles alternate such that, if betastep=1, they are
% are (0, 180) for each gamma, then (1, 179) for each gamma, up to 89, 91.
% that is, the same projected position is described twice.

%% Make the templates

disp(['Making ' num2str(nAngles) ' templates']);

tic
% allTemplates=rsMakeTemplatesQuick(angleList,map);
allTemplates=rsMakeTemplates(reshape(angleList,nAngles,3),map);
toc

% angleList=reshape(angleList,nHemiAngles,2,nGamma,3);
allTemplates=reshape(allTemplates,nt,nt,nHemi,nGamma,nHemiAngles);


%% Filter the templates according to the CTF
[nt nt nHemi nGamma nHemiAngles]=size(allTemplates);
nAngles=nHemi*nGamma*nHemiAngles;
ne=NextNiceNumber(nt*1.2);  % increase the size to allow CTF rings
ctf=meGetEffectiveCTF(mi,ne,ds);  % put in dqe, pw filter.
% evaluate a generic inverse filter
f=RadiusNorm(ne)/(mi.pixA*ds);  % frequencies for evaluating the CTF
hinv=a1./(1+(f1./f).^2)+(1-a1)*1./(1+(f0./f).^ex);
hinv(f==0)=0;
% Pad the templates to avoid ctf artifacts
% eigenSet=rsMakeEigenrefs(allTemplates,nterms,ctf.*hinv);
xTemplates=Crop(reshape(allTemplates,nt,nt,nAngles),ne,1);
nim=size(xTemplates,3);
H=ifftshift(hinv.*ctf);
msk=fuzzymask(ne,2,0.45*ne,.1*ne);
for i=1:nim
    xTemplates(:,:,i)=real(ifftn(fftn(xTemplates(:,:,i)).*H)).*msk;
end;
xTemplates=reshape(xTemplates,ne,ne,nHemi,nGamma,nHemiAngles);
ntstr=num2str(nt);
nestr=num2str(ne);
disp(['Templates expanded from ' ntstr 'x' ntstr ' to ' nestr 'x' nestr]);

%% Make the eigenreferences
eigenSet=rsMakeEigenrefs(xTemplates,nterms);
figure(2);
SetGrayscale;
plot(1-eigenSet.termVar);
% figure(1);
% ImagicDisplay(eigenSet.imgs,2);


if showTemplates
    %%  % Show the templates and the reconstructions
    timgs=reshape(eigenSet.imgs,ne*ne,nterms);
    nG=1;
    nH=2;
    nAngs=nHemiAngles;
    rImg=single(zeros(ne,ne,2*nH,nG,nAngs));
    for k=1:nHemiAngles
        for j=1:nG
            for i=1:nH
                %             rImg(:,:,2*i-1,j,k)=reshape(timgs*eigenSet.vList(:,i,j,k),ne,ne).*squeeze(eigenSet.ampList(i,j,k));
                rImg(:,:,2*i-1,j,k)=reshape(timgs*eigenSet.vList(:,i,j,k),ne,ne);
                rImg(:,:,2*i,j,k)=xTemplates(:,:,i,j,k);
            end;
        end;
    end;
    figure(1);
    ImagicDisplay2(rImg,2);
    
    
    %%  % Make a complete set of reconstructions for comparison
    %   and make the average power spectrum
    nTotal=nAngs*nGamma*nHemi;
    timgs=reshape(eigenSet.imgs,ne*ne,nterms);
    rImg=single(zeros(ne,ne,nTotal));
    vL=reshape(eigenSet.vList,nterms,nTotal);
    sp=zeros(ne,ne);
    sp0=zeros(ne,ne);
    for i=1:nTotal
        %     img=reshape(timgs*vL(:,i),ne,ne)*eigenSet.ampList(i);
        img=reshape(timgs*vL(:,i),ne,ne);
        sp0=sp0+abs(fftn(xTemplates(:,:,i))).^2;
        sp=sp+abs(fftn(img)).^2;
        rImg(:,:,i)=img;
    end;
    spr=Radial(fftshift(sp/nTotal));
    spr0=Radial(fftshift(sp0/nTotal));
    figure(2);
    semilogy([spr0 spr]);
    figure(3);
    plot([cumsum(spr0) cumsum(spr)]);
    ImagicDisplay2(rImg);
end;

% mic=rsSortVesicles(mi);  % make a copy with the vesicles sorted by position.
mic=mi;

%%  Evaluate the cc for each vesicle
disp('Evaluating cross-correlations');
figure(1); clf; SetGrayscale;
imac(imscale(GaussHP(GaussFilt(origImg,.1),.005),256,.003));
namp=1.5;
if simulateImage
    m2=m+namp*randn(n);
else
    m2=m;
end;



% mic.vesicle.ok=(mic.vesicle.s>minS) & (mic.vesicle.s<maxS);
% % mic.vesicle.s=mic.vesicle.s*0+.05;  % make the amplitudes equal.
vesInds=find(mic.vesicle.ok);

disp(['Total vesicles: ' num2str(numel(vesInds))]);
mxCC=single(zeros(n));
mxTemplInds=uint16(zeros(n));
mxVesInds=uint16(zeros(n));
mxVars=single(zeros(n));
figure(3);
SetGrayscale;
blankRadius=16;

nl=48;
msklRadius=17;  % Radius for variance averaging
mskl=fuzzymask(nl,2,msklRadius,2);  % local mask
npts=sum(mskl(:));

partRadius=10;
nves=numel(vesInds);
maskPadding=1;
for j=1:nves
    i=vesInds(j);
    maskRadii=mic.vesicle.r(i)/ds+membraneOffset*[-1 1 0]+maskPadding;
    maskRadii(3)=maskRadii(3)+abs(membraneOffset)+msklRadius+maskPadding;
    % Get the single-vesicle correlation function
    [mxVals mxInds filtImg mxValsV]=rsVesicleCorrelation3(-m2,mic,i,...
        membraneOffset,maskRadii,angleInds,eigenSet);
    ctr=round([mic.vesicle.x(i) mic.vesicle.y(i)]/ds+1);  % shift to 1-based coords
    nv=size(mxVals,1);
    
    % Do a local averaging of the squared CC in the whole vesicle.
    var=mxVals.^2;
    h=ifftshift(Crop(mskl,nv));  % average over the local mask
    filtVar=real(ifftn(fftn(var).*fftn(h))).*fuzzymask(nv,2,max(maskRadii(1:2)),.5);
    
    %     Pad to the full-sized image
    xVals=ExtractImage(mxVals,ctr,n,1);
    xTemplInds=ExtractImage(mxInds,ctr,n,1);
    xVars=ExtractImage(filtVar,ctr,n,1);
    
    %       Incorporate into the composite images
    q=xVals>mxCC;  % find maximum cc values
    mxCC(q)=xVals(q);
    mxVesInds(q)=i;  % the related vesicle indices
    mxTemplInds(q)=xTemplInds(q);  % the template indices
    
    q1=xVars>mxVars;  % find maximum variance values
    mxVars(q1)=xVars(q1);
    
    if mod(i,10)==0
        figure(3);
        imacs(mxCC);  % Show the CC function
        %         imacs(mxVars);  % Show the local variances
        title(i);
        drawnow;
    end;
end;  % loop over vesicles

partRadius=18;
save([basePath baseName 'rscc.mat'],'mxCC','mxVars','mxVesInds','mxTemplInds','m0','m1','partRadius','ds');


%% Particle finding
% This could be run as a separate program.
%
% We search for peaks in mxCC and check mxVars for excessive variance.
% MxVesInds gives the vesicle number for each particle; mxTempInds can be
% used to get the template number in case we want.
basePath='/Users/fred/matlabWork/Yunhui2/VesicleFinderData/120711/';
baseName='004_sq02_1_04';
load([basePath baseName 'rscc.mat']);
load([basePath 'Info/' baseName 'mi.mat']);
% mi=rsSortVesicles(mi);  % Sort according to position

% % partRadius=12;
% partRadius=12;
minAmp=0.75;
maxAmp=1.1;
maxVar=40;

n=size(mxCC,1);
nb=2*partRadius+4;  % size of the blanking mask
nb2=nb/2;
blankMask=1-fuzzymask(nb,2,partRadius,1);

results=zeros(1000,5);
ctrs=zeros(1000,2);
vesIndex=single(zeros(1000,1));
k=0;
nFound=0;

disp('Finding particles');
nctr=ceil((n+1/2));
mxCC2=mxCC;  % Search for cc peaks
[amp ix iy]=max2d(mxCC2);
while amp>minAmp
    ccVar=mxVars(ix,iy);
    iv=mxVesInds(ix,iy);
    result=[ iv amp round(ccVar) [ix iy]+round(nctr)];
    k=k+1;
    results(k,:)=result;
    if (amp > minAmp && amp < maxAmp && ccVar < maxVar)
        nFound=nFound+1;
        ctrs(nFound,:)=[ix iy];
        vesIndex(nFound)=iv;
        amps(nFound)=amp;
    end;
    % Blank the vicinity of the found peak
    mxCC2=Mask(mxCC2,[ix iy],blankMask);
    [amp ix iy]=max2d(mxCC2);
    if mod(k,10)==0
        imacs(mxCC2);
        title(amp);
        drawnow;
    end;
end;
results=results(1:nFound,:);
ctrs=ctrs(1:nFound,:);
for i=1:nFound
    mi.particle.x(i)=(ctrs(i,1)-1)*ds;  % should be interpolated!
    mi.particle.y(i)=(ctrs(i,2)-1)*ds;
    mi.particle.vesicle=vesIndex(1:nFound);
end;



%
figure(3);
ShowImageAndBoxes(mxCC,ctrs,20,2,[1 1 0]);

figure(2); clf;
ShowImageAndBoxes(m1,ctrs,20,2,[1 1 0],1e-3);

figure(1);
ShowImageAndBoxes(m0,ctrs,20,2,[1 1 0],1e-3);

