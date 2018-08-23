% TestMAlign
doAFlip=0;
figure(1);
SetGrayscale;
load testImgs.mat  % loads 7 images
testImgs=Crop(testImgs,48,1);
[n ny nrefs]=size(testImgs);
refs=testImgs;
if doAFlip
    refs(:,:,2)=circshift(flipud(testImgs(:,:,1)),[1 0]).*msk;  % flipped reference
end;

sigma=2;
thetas=[0 45 180 355];
thetas=0:359;
nthetas=numel(thetas);
nim=nrefs*nthetas;
imgs=single(zeros(n,n,nthetas,nrefs));
% sx=[0 0 0 0];
% sy=[0 0 0 0];

disp('Simulating images');
for i=1:nrefs
    imgs(:,:,:,i)=rsRotateImage(refs(:,:,i),thetas);
    for j=1:nthetas
        sx=j;
        sy=j/2;
       imgs(:,:,j,i)=shiftf(imgs(:,:,j,i)+sigma*randn(n,n),[-sx -sy]);  % cw rotation
    end;
end;

nImages=nrefs*nthetas

% imovie(imgs,.1);
disp('Aligning');
mode=3;
tic
xytfr=MRAlign(imgs,refs,mode,0:20:359,1,0,1);
disp('done.');
xytfr
toc
disp(' ');
% 
% [sp aliImgs]=spAlignRot(imgs,sp,refs(:,:,1));
% vals=[sp.rot sp.flip]
% 
% % [sp aliImgs]=spTransformStack(imgs, sp);
% figure(1);
% imovie(aliImgs,.1);