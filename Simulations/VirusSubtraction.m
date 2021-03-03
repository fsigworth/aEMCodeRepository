% VirusSubtraction.m
n=[720 512];
infoPath='Info/';
refMiName='Info_1/slot2_2_344_001_noDWmi.txt';
refMi=ReadMiFile(refMiName);

% cd( '/Users/fred/EMWork/Xiong/drive-download-20210209T224433Z-001/relion/work2/');

miNames=f2FindInfoFiles;
for i=1:numel(miNames)
    name=miNames{i};
    [pa,nm,ex]=fileparts(name);
mi=ReadMiFile(name);
rmi=ReadMiFile([refPath nm ex]);
m=meLoadNormalizedImage(mi,n,'m');
subplot(221);
imags(m);
drawnow;

origVesicleModel=refMi.vesicleModel;
vmSize=numel(origVesicleModel);

xvmSize=4*vmSize+1;
expVesicleModel=Crop(origVesicleModel,xvmSize);
ctrDensity=fuzzymask(xvmSize,1,xvmSize*0.4,6,1);
amp= 1;
newVesicleModel=expVesicleModel+amp*ctrDensity;
mi1=mi;
% mi1.pixA=rmi.pixA*2;
% mi1.ctf.defocus=rmi.ctf.defocus;
% mi1.ctf.deltadef=rmi.ctf.deltadef;
mi1.ctf.alpha=.02;
mi1.vesicleModel=newVesicleModel;
subplot(223);
plot(newVesicleModel);
v=meMakeModelVesicles(mi1,n);
subplot(224);
imags(v);
subplot(222);
imags(GaussFilt(m-v,.2));
pause(0.1);
WriteMiFile(mi1,name);
end;