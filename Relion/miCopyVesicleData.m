% miCopyVesicleData.m
% Copy various fields from mi files in a reference directory into the fields of
% the mi files in the working directory.

refDir='Info1/';
workDir='Info/';

d=dir(workDir);

for i=1:numel(d)
    if ~d(i).isdir && strndcmp(d(i).name,'mi.txt',6)
        workName=[workDir d(i).name];
        mi=ReadMiFile(workName);
        refName=[refDir d(i).name];
        if exist(refName,'file')
            mi1=ReadMiFile(refName);
            mi.imageSize=mi1.imageSize;
                        mi.vesicleModel=mi1.vesicleModel;
            mi.vesicle=mi1.vesicle;
            mi.log=mi1.log;
            mi.imageNormScale=mi1.imageNormScale;
            mi.imageMedian=mi1.imageMedian;
            mi.mask=mi1.mask;
            mi.timestamp=mi1.timestamp;
            WriteMiFile(mi,workName);
        else
            disp([num2str(i) '  ' Refname ' not found']);
        end;
    end;
end;
