% SimVesicleProcessing.m
% Create metadata and files for simulated micrographs containing vesicles.
% Create mi.txt metadata files and processing
% directories. The mi files can then be loaded into Vesicle_finding_GUI for
% vesicle finding, and rsRefineVesicleFitsA for refinement.

cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Data for Chris/Chris sims');
CheckAndMakeDir('Info/',1);

d=dir('Micrographs');
for i=3:numel(d)
    name=['Micrographs/' d(i).name];
    [pa,nm,ex]=fileparts(name);
    m=imread(name);
    disp(name);
    mrcName=['Merged_sm/' nm 'ms.mrc'];
    WriteMRC(m,4.24,mrcName);
    mi=MakeMiFile(name);
    mrcName1=['Merged/' nm 'm.mrc'];
    m1=Downsample(m,mi.padImageSize);
    WriteMRC(m1,mi.pixA,mrcName1);

    WriteMiFile(mi);
    mysubplot(221);
    imags(m);
    mf=GaussFilt(GaussHP(m,.01),.1);
    mysubplot(222);
    imags(mf);
    hold on;
    plot([1 1400],[450 450],'linewidth',2)
    hold off;
    mysubplot(224);
    plot(mf(:,950));
    drawnow;
end;

function mi=MakeMiFile(name);
mi0=meCreateMicrographInfoStruct14;
[imagePath,baseName,imageExtension]=fileparts(name);

imagePath=AddSlash(imagePath);

mi.baseFilename=baseName;
mi.originalBasePath=AddSlash(pwd);
mi.basePath='';
mi.moviePath='';
mi.infoPath='Info/';
mi.imageFilenames{1}=[baseName imageExtension];
mi.imagePath=imagePath;

    mi.procPath='Merged/';
    mi.procPath_sm='Merged_sm/';
mi.imageSize=4*[1440 1024];
mi.padImageSize=mi.imageSize;

    mi.pixA=1.06;

    mi.doses=50;
    mi.kV=300;
mi.camera=7;
mi.cpe=.8;
mi.weights=1;
mi.frameSets=1;

%     CTF parameters
mi.ctf=struct;
mi.ctf.defocus=2;
mi.ctf.deltadef=0;
mi.ctf.theta=0;
mi.ctf.Cs=2.7;
mi.ctf.alpha=.1;
mi.ctf.B=100;
mi.ctf.lambda=EWavelength(mi.kV);

mi.mergeMatrix=eye(3);
% Use the Grant&Grigorieff damage model
if mi.kV>250
    mi.damageModelCode='.245*f.^-1.665+2.81; % Grant&Grigorieff 300 kV';
else % 200 kV code
    mi.damageModelCode='.184*f.^-1.665+2.1; % Grant&Grigorieff 200 kV';
end;
mi.vesicle.x=[];
mi.vesicle.r=[];
mi.vesicle.s=[];
mi.vesicel.ok=[];
mi.particle.picks=[];

mi.vesicleModel=[];
mi.log=cell(0);
mi.identifier=rand;
end
