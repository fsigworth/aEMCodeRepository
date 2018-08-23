% StackClassifier.m

% Get the stack info and the stack data
[siName, pa]=uigetfile('*si.mat','Select an si file');
if isnumeric(siName)  % Cancel
    return
end;

cd(pa);
p=regexp(siName,'si.mat');
stackName=[siName(1:p-1) 'stack.mrc'];
stack=ReadEMFile(stackName);
stack=Crop(stack,48,1); %%%

[n n1 nim]=size(stack);

%% compute the mean
% Force annulus to zero
stn=NormalizeImages(stack,0,1);  % no variance normalization
stMean=mean(stn,3);
ctr=CenterOfMass(stMean);

%% Classify without alignment

mskRadius=0.35*n;
rmask=fuzzymask(n,2,mskRadius,mskRadius/10);

subplot(2,2,1);
imacs(rmask);
title('MSA mask');
drawnow;
%%
%
% Initial clssification
nfactors=20;
nclasses=50;
figure(4);
subplot(2,2,4);
[means1 nm inds fvecs]=Classifier(stn,nclasses,nfactors,2,rmask);  % 20 factors, decimate by 2.

figure(1);
ImagicDisplay(means1);
return





%%
% Initial alignment to the mean
disp('First alignment');
[xytfr, ma]=MRAlign(stn,stMean,3,2); % coarse rotation: 2 pixels at edge.
subplot(2,2,3);
plot(xytfr(1),xytfr(2));
title(shifts);
subplot(2,2,2);
q=mean(ma,3);
imacs(q);
title('Totsum');

mskRadius=0.35*n;
rmask=fuzzymask(n,2,mskRadius,mskRadius/10);

subplot(2,2,1);
imacs(rmask);
title('MSA mask');
drawnow;
return
%%
%
% Initial clssification
nfactors=20;
nclasses=100;
figure(4);
subplot(2,2,4);
[means1 nm inds fvecs]=Classifier(ma,nclasses,nfactors,2,rmask);  % 20 factors, decimate by 2.

figure(1);
ImagicDisplay(means1);

%%
%  -------- Ask the user for reference numbers ---------
refpicks=input('Reference picks, e.g. [23 26]? ');
if numel(refpicks)<1
    refpicks=[23 26];
end;
refpicks

refis=means1(:,:,refpicks);
tstep=2;  % quantum of rotation

disp('Second alignment');
transSD=10;
alprior=-Radius(n).^2/(2*transSD^2);
[xytfr ali]=MRAlign(ma,refis,3,tstep,alprior);

% Now that the images are aligned, do a second classification
nclasses=100;
nfactors=30;
disp('Sharp filter');
ali=SharpFilt2(ali,.45,.05);

amask=fuzzyblob(64,[20 0 -2],4);  % elliptical mask
amask=rmask; %%%%--------use simple mask instead.-----------
figure(4);
subplot(2,2,1);
imacs(amask);
title('Classification mask');
drawnow;

figure(2);
mskmeans=means1.*repmat(amask,[1,1,size(means,3)]);
ImagicDisplay(mskmeans);
drawnow;


disp('Classifier');
[means nm inds fvecs eigims sval]=Classifier(ali,nclasses,nfactors,1,amask);  % 20 factors, decimate by 2.

figure(6);
ImagicDisplay1(means,1);
figure(2);
ImagicDisplay1(eigims,1);
% figure(5);
% bar(sval);
%%
outPath='';

WriteImagic(means,[outPath 'ClassMeans']);
WriteImagic(eigims,[outPath 'EigenImages2']);
WriteImagic(refis,[outPath '2Refs']);

