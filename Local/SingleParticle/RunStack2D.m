% RunStack2D
% script to do 2D classification
clear;  % clear out all variables.
nEigs=15;
nClasses=20;

stackPath='/Volumes/TetraData/EMWork/Hideki/121118/PC1PC2liposome_SUP_14k/Merged/';
stackName='sq04_20024stack64.mat';

load([stackPath stackName]);
% loads stack, ctfs and si
si.pixA=5.8;
n=size(stack,1);

% then center the particles.
disp('NormalizeAndCenter');
[sp st1]=spNormalizeAndCenter(stack,si.pixA,.48,2);
limL=mean(sp.amp)-2*std(sp.amp);
limU=mean(sp.amp)+2*std(sp.amp);
limL=4;%%%%
st1s=st1(:,:,sp.amp>limL & sp.amp<limU);
%%
disp('Invert and filtering');
st1=-GaussFilt(st1s,.3,1);
disp('MakeEigenimages');
msk=fuzzymask(n,2,n*.3,n*.1);
eigs=spMakeEigenimages(st1,nEigs,1,msk);

imats(eigs,2)  % show the eigenimages scaled up by 2
%%
disp('Classify');
[sp means1]=spClassify(st1,sp,eigs,nClasses);  % make 50 class means
imats(means1,2)
%%
ref=means1(:,:,13);
% ref=fuzzymask(32,2,[10 5],4);  % artificial elliptical reference
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
imats(means2,2)
%%
%  pick some of the class means as new references
refs=means2(:,:,[7 11]);
imats(refs,2);
for i=1:size(refs,3)
    refs(:,:,i)=refs(:,:,i).*msk;
end;


disp('Second alignment');
[sp aliImgs2]=spAlignRotTrans(st1,sp,refs);
%%
eigs=spMakeEigenimages(aliImgs2,nEigs*2,1,msk);
[sp means3]=spClassify(aliImgs2,sp,eigs,nClasses);
imats(means3,2);
