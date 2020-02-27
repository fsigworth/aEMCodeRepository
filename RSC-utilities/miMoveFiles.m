% miMoveFiles.m
% Move all the mi.txt and corresponding mie.mat files to another directory
doExec=1;
infoDir='Info_small2/';
outDir='Info_small22/';

miNames=f2FindInfoFiles(infoDir);
nn=numel(miNames);
if nn<1
    return
end;

CheckAndMakeDir(outDir,1);

for i=1:nn
    miName=miNames{i};
    [pa,nm,ex]=fileparts(miName);
    pax=AddSlash(pa);
    str1=['mv ' miName ' ' outDir nm ex];
    mieNm=[nm 'e.mat'];
    if exist([pax mieNm],'file')
        str2=['mv ' pax mieNm ' ' outDir mieNm];
    else
        str2='';
        disp(['not present: ' mieNm]);
    end;
    if doExec
        disp(str1);
        system(str1);
        disp(str2);
        system(str2);
    end;
end;  