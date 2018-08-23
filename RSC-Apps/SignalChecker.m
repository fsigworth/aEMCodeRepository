% SignalChecker
% Estimates the vesicle signal in a micrograph.
% 
defaultPixA=1.7;
initialDefocus=1;
cpe=16;
iCamera=1;
minAmp=1e-3;
maxAmp=5e-3;

% Have the user select an image files
[fname, pa]=uigetfile('*.*','Select an image file');
cd(pa);
fnames{1}=fname;
%%
figure(1);
clf;
SetGrayscale;


mi=meCreateMicrographInfoStruct12;
[m, pixA, mi.doses]=meReadImagesNorm(fnames, cpe, 0,defaultPixA);
dose=mi.doses
mi.pixA=pixA;
mi.imageSize=size(m);
mi.weights=1;
mi.kV=200;

modelSpectrum=CCDModelSpectrum2D(iCamera);  % handle DE-12 or CCD
mp=mePreWhiten(m,modelSpectrum);
mi.ctf=meFitCTF(mp,mi.pixA,initialDefocus);  % Actual CTF fitting
defocus=[mi.ctf.defocus mi.ctf.deltadef]

vLipid=1.6;
thk=60;
rise=6;
% Create the model, which is sampled in units of the original pixel size.
nm0=ceil(30/pixA)*2+1;  % array for vesicle model; 60A nominal
mi.vesicleModel=fuzzymask(nm0,1,thk/pixA/2,rise/pixA)...
    *vLipid;  % units of V
%%
rPars=[100 300 10];
meFindVesicles3('end');  % deallocate the finder.
[mi1, t]=meFindVesicles3(mp, mi, rPars);
%%
minAmp=1e-3;
maxAmp=5e-3;
vesicleAmps=[minAmp maxAmp];
mins=inf;
nVesOld=0;
% Loop through finding groups of 50 vesicles.
while mins>minAmp
    [mi1, t]=meFindVesicles3('next',50,vesicleAmps);
    subplot(1,2,1);
    imacs(t.ms-t.umodel);
    mins=t.globalmax;
    nves=numel(mi1.vesicle.s);
    %     Make the scatterplot of vesicle radii and amplitudes
    subplot(1,2,2);
    plot(mi1.vesicle.r(1:nves)*mi1.pixA,1000*mi1.vesicle.s(1:nves),'b.');
    xlabel('Vesicle radius, Å');
    ylabel('Amplitude');
    title(nves);
    drawnow;
    %     Exit the loop when we can't find any more vesicles
    if nves<=nVesOld
        break
    end;
    nVesOld=nves;
end;
disp('done.');
