% CopyMiVesicles.m
% Copy mi.vesicle field from info files in one directory to another.

cd '~/project/20181025/';
rootDir=pwd;
inDir='80Frames/';
outDir='20Frames/';
frameSets=[2 21];
cd(inDir);
names=f2FindInfoFiles;
nmi=numel(names);
cd (rootDir);

for i=1:nmi
    mi1=ReadMiFile([inDir names{i}]);
    outFile=[outDir names{i}];
    mi2=ReadMiFile(outFile);
    mi2.vesicleModel=mi1.vesicleModel;
    mi2.vesicle=mi1.vesicle;
    % Correct the dose information
    mi2.frameSets=frameSets;
    mi2.doses=sum(mi1.frameDose(frameSets(1,1):frameSets(1,2)));

    WriteMiFile(mi2,outFile);
    disp(outFile);
end;

