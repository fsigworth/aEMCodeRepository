% meSaveSingleImages
% Given an info structure mi and the original images, write out phase-
% flipped individual images and separate mi files.  I use this to try
% refining vesicle fits to each separate image.
% fs 16 feb 13
mcDS=2;
cpe=15;
overridePixA=2.9;
iCamera=1;  % F20 US4000
modelSpectrum=CCDModelSpectrum2D(iCamera);  % handle DE-12 or CCD

[fname pa]=uigetfile('*mi.mat','Select an mi file');
[rootPath infoPath]=ParsePath(pa);
cd(rootPath);
disp(['Reading ' fname]);
load([infoPath fname]);

fnames=cell(0);
for j=1:nim
    fnames{j}=[mi.basePath mi.imagePath mi.imageFilenames{j}];
    disp(mi.imageFilenames{j});
end;
%%
% Read the images
[m mi.pixA mi.doses]=meReadImages(fnames,cpe,0,overridePixA);
sz=size(m);
mi.imageSize=sz(1:2);
% operate with the pre-whitening filter
%     disp('mePreWhiten');
m=mePreWhiten(m,modelSpectrum);

miOrig=mi;
nim=numel(miOrig.doses);
for j=1:nim
    mi=miOrig;
    mi.ctf=mi.ctf(j);
    mi.doses=mi.doses(j);
    if j==1
        mc0=m(:,:,1);
    else
        mc0=AffineTransform(m(:,:,j),miOrig.mergeMatrix(:,:,j));
    end;
    mc1=Downsample(mc0,mi.imageSize/mcDS);
    theCTF=CTF(mi.imageSize/mcDS,mi.pixA*mcDS,mi.ctf);
    %                 --doesn't include pre-dose effects!
    %                 phase-flip
    mc=-real(ifftn(fftn(mc1).*ifftshift(sign(theCTF))));
    mi.baseFilename=[mi.baseFilename 'exp' num2str(j)];
    procname=[mi.basePath mi.procPath mi.baseFilename];  % Use the processed image directory
    WriteMRC(mc,mcDS*mi.pixA,[procname 'm.mrc']);
    infoFullName=[basePath infoPath mi.baseFilename 'mi.mat'];
    save(infoFullName,'mi');
    disp(['Info file written: ' infoFullName]);
end;
