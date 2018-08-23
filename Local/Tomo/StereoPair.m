% function StereoPair
% Align a tilt pair
% Images are always displayed after rot90 (ccw rotation) to allow
% cross-eyed viewing.

% ----------------Get merged micrographs---------------

angle=15;
ds1=2;
ds2=1;

% Have the user select some info files
[gMiNames, pathName]=uigetfile('*mi.*','Select a pair of mi files','multiselect','on');
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
mi0=ReadMiFile([gInfoPath gMiNames{ind0}]);
mi=ReadMiFile([gInfoPath gMiNames{ind1}]);
%         if (isfield(mi,'tiltAngles') && numel(mi.tiltAngles)>1 && mi.tiltAngles==0) && ~forceTransform
%             disp(['Skipping ' gMiNames{ind1}]);
%         else
mi0.basePath=rootPath;

disp('Reading untilted, merged image:');
[m0, mergedImgName]=meReadMergedImage(mi0,1);
disp(['    ' mergedImgName]);
ds0=mi0.imageSize(1)/size(m0,1);
%
% For the tilted image, re-merge it
disp('Reading tilt image:');
mi.basePath=rootPath;
[mRaw, mi.pixA, mi.doses]=meReadImagesNorm(mi,mi.cpe,0,mi.pixA);
% [m1,tiltImgName]=meReadMergedImage(mi,1);
disp(['     ' mi.baseFilename]);

disp('Prewhitening');
modelSpectrum=CCDModelSpectrum2D(mi.camera);  % handle DE-12 or CCD
mRaw=mePreWhiten(mRaw,modelSpectrum);
    [~, m1]=meCombineImages2(mRaw,mi,ds0,1);
    m1=-m1;  % go back to reversed contrast
    mskSize=round(mi.imageSize(1)/4);
    msk=meMakeMergedImageMask(mskSize,mi.mergeMatrix);
    mi=meInsertMask(msk,mi,1);  % merge mask is the first one in the mask stack.

pixA=mi.pixA*ds0;

%
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
p=[95 angle 1];  % Tilt axis, tilt, magnification
pSteps=[1 2 0.02];
pSteps2=pSteps/20;
shift=[0 0];
niters=80;
niters=35;
hpIter=40;  % No. of iterations for HP filter to get turned up to max
prior=(1e-5*RadiusNorm(nw).^2);
rs.mask=fuzzymask(nw,2,nw*.4,nw/8);
im0=stbin(:,:,i0);
me0=mean(im0(:));  % Untilted image
im0=im0-me0;

msk=SquareWindow(nw,nw/16);
imt=stbin(:,:,it).*msk;
rs.fc=.08;  % maximum hp filter, scaled by ds
rs.f0=.2;  % filtering for display
rs.ft=.1;  % for display also
rs.ds=ds1*ds2;
rs.pixA0=pixA;

p=Simplex('init',p,pSteps,[1 1 1]);
%         shStack=stbin;
        rs.mxShift=1;
        rs.prior=prior/10;
        rs.lpf=.05;  % lowpass of cc peak, not scaled by ds
for i=1:min(hpIter,niters)
    %%
    if i==8
        rs.prior=prior;
        rs.mask=fuzzymask(nw,2,nw*.2,nw/8);

    elseif i==20
        rs.mxShift=.4;
         rs.lpf=.15;
    end;
    rs.hpf=min(i/hpIter*rs.fc,rs.fc);
    rs.mxShift=0.8+(i<10)*.5;
    [p,shift]=StereoFitIteration(p,shift,im0,imt,rs);
    title(i);
end;
%  p=Simplex('init',p,pSteps,[1 1 1]);
for i=hpIter+1:niters
     [p,shift]=StereoFitIteration(p,shift,im0,imt,rs);
      title(i);
end;   
%         Make the final matrix
p=Simplex('centroid');
disp([p shift]);
[~,~,T,imtc]=StereoFitIteration(p,shift,im0,imt,rs);
subplot(224);
imacs(rot90(im0-imtc/mi.doses(1)));
drawnow;

rs.mask=fuzzymask(n0,2,n0/4,n0/8);

mi.mergeMatrix=T;
mi.tiltAxis=p(1);
mi.tiltAngles=p(2);
fullName=meSaveMiFile(mi);
disp(['mi file saved: ' fullName]);

%%  Do the nice transform on the full-sized image and save the result
% [~, m10]=meCombineImages2(mRaw,mi,1,1);
% mc10=AffineTransform(-m10,T);
mc10=AffineTransform(m1,T);
mc=Downsample(mc10,size(m0));
procname=[mi.basePath mi.procPath mi.baseFilename];  % Use the processed image directory
jpegname=[mi.basePath mi.procPath 'jpeg/' mi.baseFilename];

% For jpeg files we trim .1% of values from the intensity
% histograms.
WriteMRC(mc,pixA,[procname 'tm.mrc']);
disp(['Wrote transformed image: ' procname 'tm.mrc']);
WriteJpeg(mc,[jpegname 'tm.jpg']);
disp(['Wrote transformed image: ' jpegname 'tm.jpg']);
% imacs(rot90(mc));

%         end;  % if skipping
%     end; % for j
% end; % for iHole
%
