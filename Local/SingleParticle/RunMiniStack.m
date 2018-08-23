% RunMiniStack
% script to do 2D classification as demonstrated in class 4/10
clear;  % clear out all variables.

% First, you have to select File -> Import data,
%  which loads the "ministack.dat" file

whos % check to see that rawStack has been loaded

% then center the particles.
[sp st1]=spNormalizeAndCenter(rawStack,5.8);
%%
eigs=spMakeEigenimages(st1,20);

imats(eigs,2)  % show the eigenimages scaled up by 2
%%
[sp means1]=spClassify(st1,sp,eigs,50);  % make 50 class means
imats(means1,2)
%%
ref=means1(:,:,8);
% ref=fuzzymask(32,2,[10 5],4);  % artificial elliptical reference
clf
imacs(ref)
%%
[sp aliImgs]=spAlignRotTrans(st1,sp,ref);
imacs(mean(aliImgs,3):  % show the mean of all the aligned images.
%%
eigs2=spMakeEigenimages(aliImgs,20);
[sp means2]=spClassify(aliImgs,sp,eigs,50);
imats(means2,2)
%%
%  pick some of the class means as new references
refs=means2(:,:,[1 50 11]);
imats(refs,2)
%%
[sp aliImgs2]=spAlignRotTrans(st1,sp,refs);
%%
eigs=spMakeEigenimages(aliImgs2,20);
[sp means3]=spClassify(aliImgs2,sp,eigs,100);
imats(means3,2);
