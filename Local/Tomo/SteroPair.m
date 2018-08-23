% function StereoPair
% Align a tilt pair

% ----------------Get merged micrographs---------------

angle=15;
ds1=2;
ds2=1;

% Have the user select some info files
[gMiNames, pathName]=uigetfile('*mi.mat','Select a pair of mi files','multiselect','on');
if isnumeric(pathName) % File selection was cancelled
    return
end;
if ~iscell(gMiNames)
    gMiNames={gMiNames};
end;
[rootPath, gInfoPath]=ParsePath(pathName);
cd(rootPath);


%%
% nHoles=ceil(numel(gMiNames)/14);
% nHoles=1;
% for iHole=1:nHoles
%     nShots=min(7,numel(gMiNames)-14*(iHole-1)-7);
%     for imgIndex=1:nShots
ind0=1;  % zero-tilt image
ind1=2;  % tilted image
mi0=load([gInfoPath gMiNames{ind0}]);
mi0=mi0.mi;
mit=load([gInfoPath gMiNames{ind1}]);
mi=mit.mi;
%         if (isfield(mi,'tiltAngles') && numel(mi.tiltAngles)>1 && mi.tiltAngles==0) && ~forceTransform
%             disp(['Skipping ' gMiNames{ind1}]);
%         else
mi0.basePath=rootPath;

[m0, mergedImgName]=meReadMergedImage(mi0,1);
disp('Reading untilted, merged image:');
disp(mergedImgName);
ds0=mi0.imageSize(1)/size(m0,1);
%%
% For the tilted image, re-merge it
disp('Reading tilt image:');
[mRaw, mi.pixA, mi.doses]=meReadImagesNorm(mi,mi.cpe,0,mi.pixA);
disp('Prewhitening');
modelSpectrum=CCDModelSpectrum2D(mi.camera);  % handle DE-12 or CCD
mRaw=mePreWhiten(mRaw,modelSpectrum);
    [~, m1]=meCombineImages2(mRaw,mi,ds0,1);
    m1=-m1;  % go back to reversed contrast
    mskSize=round(mi.imageSize(1)/4);
    msk=meMakeMergedImageMask(mskSize,mi.mergeMatrix);
    mi=meInsertMask(msk,mi,1);  % merge mask is the first one in the mask stack.

mi.basePath=rootPath;
pixA=mi.pixA*ds0;

st=single(m0);
st(:,:,2)=m1;
[nx, ny, nim]=size(st);
n0=[nx ny];
%
figure(1);
SetGrayscale;
disp('Binning');
stbin=zeros([n0/ds1 nim],'single');
nw=round(floor(n0/ds1)/ds2);
for i=1:nim
    q=single(st(:,:,i));
    me=mean(q(:));
    stbin(:,:,i)=Downsample(BinImage(q-me,ds1),nw);
end;
it=2;  % tilted image index
i0=1;  % untilted image index
%         kp=2;
%         qpars=zeros(nim,4);
p=[90 angle 1];  % Tilt axis, tilt, magnification
pSteps=[5 5 0.1];
shift=[0 0];
niters=50;
hpIter=20;  % No. of iterations for HP filter to get turned up to max
rs.prior=(1e-5*RadiusNorm(nw).^2);
rs.mask=fuzzymask(nw,2,nw/4,nw/8);
im0=stbin(:,:,i0);
me0=mean(im0(:));  % Untilted image
im0=im0-me0;

imt=stbin(:,:,it);

rs.fc=.12;
rs.f0=.2;
rs.ft=.1;
rs.ds=ds1*ds2;
rs.pixA0=pixA;

p=Simplex('init',p,pSteps,[1 1 0]);
%         shStack=stbin;
for i=1:hpIter
    %%
    rs.hpf=min(i/hpIter*rs.fc,rs.fc);
    [p,shift]=StereoFitIteration(p,shift,im0,imt,rs);
    title(i);
end;
p=Simplex('init',p,pSteps,[1 1 0]);
for i=hpIter+1:niters
     [p,shift]=StereoFitIteration(p,shift,im0,imt,rs);
      title(i);
end;   
%         Make the final matrix
p=Simplex('centroid');
disp([p shift]);
[~,~,T,imtc]=StereoFitIteration(p,shift,im0,imt,rs);
subplot(224);
imacs(im0-imtc/mi.doses(1));
drawnow;

rs.mask=fuzzymask(n0,2,n0/4,n0/8);

mi.mergeMatrix=T;
mi.tiltAxis=p(1);
mi.tiltAngles=p(2);
fullName=meSaveMiFile(mi);
disp(['mi file saved: ' fullName]);

%%  Do the nice transform on the full-sized image and save the result
[~, m10]=meCombineImages2(mRaw,mi,1,1);
mc10=AffineTransform(-m10,T);
mc=Downsample(mc10,size(m0));
procname=[mi.basePath mi.procPath mi.baseFilename];  % Use the processed image directory
jpegname=[mi.basePath mi.procPath 'jpeg/' mi.baseFilename];

% For jpeg files we trim .1% of values from the intensity
% histograms.
WriteMRC(mc,pixA,[procname 'm.mrc']);
disp(['Wrote transformed image: ' procname 'm.mrc']);
WriteJpeg(mc,[jpegname 'm.jpg']);
disp(['Wrote transformed image: ' jpegname 'm.mrc']);
imacs(mc);

%         end;  % if skipping
%     end; % for j
% end; % for iHole
%
