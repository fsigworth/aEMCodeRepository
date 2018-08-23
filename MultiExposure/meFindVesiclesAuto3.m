% meFindVesiclesAuto3
% Test for meFindesicles3.m
% 
% Automatic vesicle picker.
% Asks for the location of an mi file, assumed to be in
%   --basePath-- Info/    (assuming mi.infoPath = 'Info/')
%  and then finds a merged image in
%   --basePath-- Merge/   (assuming that mi.procPath = 'Merge/')


fileSuffix1='mc.mrc';
fileSuffix1='m.jpg';

%     [fname inPath]=uigetfile('*mi.mat','Select an mi file');
fname='004_sq02_1_04mi.mat';
inPath='/Volumes/TetraData/EMWork/Hideki/120711/AMPA_R_dialyzed_centrifuged_sampleA/Info/';
cd(inPath)
load(fname);
basePath='../';
% end;

%%

iname1=[basePath mi.procPath mi.baseFilename fileSuffix1];
iname2=[basePath mi.procPath mi.baseFilename 'm.mrc'];
if FileExists(iname1)
    m=single(ReadEMFile(iname1));
    iname=iname1;
elseif FileExists(iname2)
    m=single(ReadEMFile(iname2));
    iname=iname2;
else
    error(['Can''t find merged image in this directory: ' pwd mi.procPath]);
end;

% membrane model
vLipid=1.6;
thk=60;
rise=6;
% Create the model, which is sampled in units of the original pixel size.
nm0=ceil(30/mi.pixA)*2+1;  % array for vesicle model; 60A nominal
mi.vesicleModel=fuzzymask(nm0,1,thk/mi.pixA/2,rise/mi.pixA)...
    *vLipid*mi.pixA;  % units of V.A per voxel

figure(2);
SetGrayscale;


rPars=[50 500 10];

mi1=meFindVesicles3(m, mi, rPars);
%%
mins=1;
nVesOld=0;
while mins>.01
    [mi1 t]=meFindVesicles3('next',20,[.025 .06]);
    mins=min(mi1.vesicle.s)
    nves=numel(mi1.vesicle.s)
    imacs(t.umodel);
    if nves<=nVesOld
        break
    end;
    nVesOld=nves;
end;



return

% [mi t]=meFindVesicles3(m, mi, rPars)          --initialize vesicle finder
% mi = meFindVesicles3('next',maxN,minThresh);  --find up to maxN vesicles.
% mi = meFindVesicles3(m, mi) -- new image, old parameters
% meFindVesicles('end');      --deallocates the persistent variables.
% t is a persistent structure that contains information for the

% figure(1);
% SetGrayscale;
% subplot(2,3,1);
% imacs(ms);

return
%%

disp('Final subtraction of vesicles');
vm=meMakeModelVesicles(mi,size(m));
mv=m-vm;
%%
figure(2);
imacs(GaussFilt(mv,.2));

figure(1);
subplot(2,3,1);
imacs(BinImage(m,4));
subplot(2,3,2);
imacs(vm);
subplot(2,3,3);
plot(mi.vesicle.r*mi.pixA,mi.vesicle.s,'k.','markersize',10);

subplot(2,3,4);
imacs(BinImage(mv,4));
title('Subtracted');
subplot(2,3,5);
hist(mi.vesicle.s,50);
xlabel('Vesicle amplitude s');
drawnow;

mi.basePath=ParsePath(inPath);  % make it the local path
bname=[mi.basePath mi.procPath mi.baseFilename];
WriteMRC(mv,ds0*mi.pixA,[bname 'mv.mrc']);
jname=[mi.basePath mi.procPath 'jpeg/' mi.baseFilename];
imwrite(uint8(imscale(rot90(mv),256,1e-3)),[jname 'mv.jpg']);

