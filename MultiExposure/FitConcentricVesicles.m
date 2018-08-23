% Fit concentric vesicles
% Try fitting the tiny interior vesicles in the PC1-PC2 prep
% fs 7 Dec 12
basePath='/Volumes/TetraData/EMWork/Hideki/121201/PC1PC2liposome_14kSUP/';
load([basePath 'Info/sq02_10000mi.mat']);
m=ReadEMFile([basePath mi.procPath mi.baseFilename 'm.mrc']);
indices=find(mi.vesicle.s>.05 & mi.vesicle.s<.1);
v=meMakeModelVesicles(mi,size(m),indices);
msub=m-v;
%%
figure(2);
SetGrayscale;
imacs(GaussFilt(msub,.1));
drawnow;
%%
ind0=max(indices)+1;
ind=ind0;
mi.vesicle.x(ind)=2930;
mi.vesicle.y(ind)=2475;
mi.vesicle.r(ind)=24;
mi.vesicle.s(ind)=.08;

ind=ind+1;
mi.vesicle.x(ind)=708;
mi.vesicle.y(ind)=2124;
mi.vesicle.r(ind)=25;
mi.vesicle.s(ind)=.08;

ind=ind+1;
mi.vesicle.x(ind)=2527;
mi.vesicle.y(ind)=2715;
mi.vesicle.r(ind)=20;
mi.vesicle.s(ind)=.08;


v2=meMakeModelVesicles(mi,size(m),ind0:ind);

filtSub=GaussFilt(msub-v2,.1);
fss=imscale(filtSub,256,1e-3);
imac(fss);

outFile=[basePath mi.procPath 'jpeg/' mi.baseFilename 'msv.jpg'];
imwrite(uint8(imscale(rot90(msub-v2),256,1e-3)),[outFile]);
