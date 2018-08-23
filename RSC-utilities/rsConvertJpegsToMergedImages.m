% rsConvertJpegsToMergedImages.m

ds=2;
procPath='Merged/';

   [names, pathName]=uigetfile({'*mi.txt';'*mi.mat'},'Select mi files','multiselect','on');
    if isnumeric(pathName) % File selection was cancelled
        return
    end;
    if ~iscell(names)
        names={names};
    end;
    [rootPath, infoPath]=ParsePath(pathName);
    cd(rootPath);
    miNames=cell(0);
    for i=1:numel(names)
        miNames{i}=[infoPath names{i}];
    end;

CheckAndMakeDir(procPath,1);

    nim=numel(miNames);
    for i=1:nim
        mi=ReadMiFile(miNames{i});
        jName=['jpeg/' mi.baseFilename 'm.jpg'];
        m=single(ReadEMFile(jName));
        
%         Scale to match merged image scale
        sp1=RadialPowerSpectrum(m);
        n2=numel(sp1);
        s1=sqrt(mean(sp1(round(n2*.75):round(n2*.95))));
        s2=sqrt(mi.doses(1)*ds*mi.pixA);
        ms=m/(s1*s2);
        outName=[procPath mi.baseFilename 'mz.tif'];
        disp(outName);
        WriteZTiff(ms,mi.pixA*ds,outName);
    end;