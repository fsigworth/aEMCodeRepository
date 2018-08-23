% meMakeMiForPPImage

[fname pa]=uigetfile('*.*','Select a micrograph');
if numel(fname)<4 % nothing selected
    return
end;
cd(pa);
fileName=fname;
disp(fileName);
%%

[m0 pixA]=ReadEMFile(fileName);
if pixA==0
    pixA=input('Enter pixel size in angstroms: ');
end;

n=size(m0);
m1=RemoveOutliers(m0);
df=1/(n(1)*pixA/10);  % frequency step in nm^-1

%%
cpe=150;

m=m1/cpe;
doses=mean(m(:));
m=m-doses;


%%
mi=meCreateMicrographInfoStruct12;
[pa nm ex]=fileparts(fileName);
mi.baseFilename=nm;
mi.basePath=AddSlash(pa);
mi.imagePath='';
mi.procPath='Merge/';
mi.infoPath='Info/';
mi.imageFilenames{1}=fileName;
mi.imageSize=n;
mi.pixA=pixA;
mi.doses=doses;

ctf.lambda=EWavelength(200);
ctf.defocus=0;
ctf.deltadef=0;
ctf.theta=0;
ctf.alpha=1;
ctf.Cs=7;
ctf.B=250;
ctf.cuton=.001;
mi.ctf=ctf;


fc=.001*pixA;
mf=SharpHP(m,fc,fc/10);

WriteMRC(mf,mi.pixA,[mi.basePath mi.procPath mi.baseFilename 'm.mrc']);
WriteJpeg(mf,[mi.basePath mi.procPath mi.baseFilename 'm.jpg']);

save([mi.basePath mi.infoPath mi.baseFilename 'mi.mat'],'mi');



