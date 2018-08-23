% TestVesicleRemoval.m

% Generic inverse filter params
f0=.02; f1=.0015; a1=.4; ex=4;
figure(1); SetGrayscale;

basePath='/Users/fred/matlabWork/Yunhui2/VesicleFinderData/120711/';
miName='Info/004_sq02_1_04smi.mat';
miName='Info/007_sq02_1_07smi.mat';
% miName='Info/010_sq02_1_10smi.mat';
load([basePath miName]);
baseName=['Merged/' mi.baseFilename];
m=ReadEMFile([basePath baseName 'm.mrc']);
m=BinImage(m,2); %%%%%
n=size(m);
ds=mi.imageSize(1)/n(1);
imacs(m);
drawnow;

if max(mi.vesicleModel)>2
    mi.vesicleModel=mi.vesicleModel/mi.pixA;
end;
disp('Making model vesicles');
tic
vim=meMakeModelVesicles(mi,n);
toc

imacs(m-vim);
title([baseName 'm.mrc'],'interpreter','none');
drawnow;
%
outName=[basePath baseName 'mv.mrc'];
% WriteMRC(m-vim,mi.pixA*ds,outName);


