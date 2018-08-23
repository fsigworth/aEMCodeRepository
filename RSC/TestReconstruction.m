% TestReconstruction

load /Volumes/TetraData/Structures/AMPAR/3KG2RotMap5.8A.mat
map=Crop(map/2000,48);  % approx normalization to give s/n=10, pad to 48 pixels.

% ri is the run info structure
n=size(map,1);
ri.n=n;
ri.membraneOffset=-12;  % this map is downsampled by 2.
ri.pixA=5.8;

ri.angleStep=10;  % degrees.  5 degrees would be more realistic
ri.symmetry=2;
ri.nGamma=round(90/ri.angleStep);
ri.maxTranslation=15;
sigmaN=1;

angleList=rsMakeTemplateAngles(ri);

alphaSteps=4;
alphas=(-alphaSteps:alphaSteps)*ri.angleStep;

% make templates
disp(['Making ' num2str(size(angleList,1)) ' templates']);
templates=rsMakeTemplates(angleList,map);
disp('done');

% make an image and a set of rotated copies.
trans=[-2 0];
img=circshift(rsMakeTemplates([0 90 0],map),trans)...
    +fuzzymask(n,2,n*0.45,n*.05)*sigmaN.*randn(n,n);

rotImgs=rsRotateImage(img,alphas);

figure(1);
ImagicDisplay(templates);
figure(2);
ImagicDisplay(rotImgs,2);

% model parameters
model.noiseSD=1;
model.clickSD=1; % pixels
model.rockSD=2;  % degrees
model.bobSD=1;   % pixels
model.pIO=0;
model.b0=0;

