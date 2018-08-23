% AnchiRunStack2D
% script to do 2D classification
% For Anchi's movie

clear;  % clear out all variables.
nEigs=10;
nClasses=10;

stackPath='/Volumes/TetraData/EMWork/Anchi/12apr26c/Ptcls/';
stackName='12apr26c008ed-ali.img';
pixA=1.8;
maskRadius=.38;

stack=ReadImagic([stackPath stackName]);

n=size(stack,1);

% then center the particles.
disp('NormalizeAndCenter');
[sp st1]=spNormalizeAndCenter(stack,pixA,.48,2);
% limL=mean(sp.amp)-2*std(sp.amp);
% limU=mean(sp.amp)+2*std(sp.amp);
% limL=4;%%%%
% % st1s=st1(:,:,sp.amp>limL & sp.amp<limU);
%%
disp('Invert and filtering');
st1=-GaussFilt(st1,.3,1);
disp('MakeEigenimages');
msk=fuzzymask(n,2,n*maskRadius,n*.1);
eigs=spMakeEigenimages(st1,nEigs,1,msk);

imats(eigs,1)  % show the eigenimages scaled up by 2
%%
disp('Classify');
[sp means1]=spClassify(st1,sp,eigs,nClasses);  % make class means
imats(means1)
%%
ref=means1(:,:,3);

clf
imacs(ref)
%%
disp('Align');
[sp aliImgs]=spAlignRotTrans(st1,sp,ref.*msk);
imacs(mean(aliImgs,3));  % show the mean of all the aligned images.
drawnow;
%%
disp('Second classification');
eigs2=spMakeEigenimages(aliImgs,nEigs);
[sp means2]=spClassify(aliImgs,sp,eigs,nClasses);
imats(means2)
%%
%  pick some of the class means as new references
refs=means2(:,:,[3 9]);
imats(refs);
for i=1:size(refs,3)
    refs(:,:,i)=refs(:,:,i).*msk;
end;


disp('Second alignment');
[sp aliImgs2]=spAlignRotTrans(st1,sp,refs);
%%
eigs=spMakeEigenimages(aliImgs2,nEigs*2,1,msk);
[sp means3]=spClassify(aliImgs2,sp,eigs,nClasses);
imats(means3);
%%

figure(2);
SetGrayscale;
s=RadialPowerSpectrum(aliImgs2);
semilogy(s)
%%
nim=size(aliImgs2,3);
spect=zeros(n,n);
for i=1:nim
    spect=spect+abs(fftn(aliImgs2(:,:,i))).^2;
end;
imacs(fftshift(spect.^.5));
