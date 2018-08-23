% FilterMergedImagesLFFlat.m
LFAmp=0.4;
if ~exist('uiPath','var')
    infoPath='/Volumes/TetraData/EMWork';
end;
cd(infoPath);
[filename infoPath]=uigetfile('*.mat','Select an info file');
infos=meLoadInfoFiles(infoPath);
nim=numel(infos);
if nim<1
    error(['No info files found in ' infoPath]);
end;
basePath=ParsePath(infoPath);

for i=1:nim
    mi=infos{i};
    mergeFilename=[basePath mi.procPath mi.baseFilename 'm.mrc'];
    disp(['Reading ' mergeFilename]);
    m=ReadEMFile(mergeFilename);
    [mf c T mfg]=meFilterImageFlatLF(m,mi,LFAmp);
    outDir=[basePath mi.procPath 'filteredL/'];
    if ~DirectoryExists(outDir)
        mkdir(outDir);
    end
    outName=[outDir mi.baseFilename 'mlf.jpg'];
    disp(['Writing ' outName]);
    WriteJpeg(mfg,outName,1e-3);
end;
