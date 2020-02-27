% miSortByJpegs

%  We assume we have allNames and allMis loaded, e.g. by running MiLoadAll,
% and we're in the expt directory.
% We create directories like Info_0010/ and put the mi files of the
% appropriate indices into them.
allNames=f2FindInfoFiles;
mPath='/Volumes/Drobo4/191228.8/Merged_redo/';
mergePath='Merged/';

%
doWrite=0;
writeMicrographs=1;
% Indices for 191228.7
% indices= ...
% [ 256 301
% 511 560
% 720 813
% 1174 1223
% 1331 1376
% 1459 1460
% 1637 1733
% 2040 2090
% 2196 2292 ];

% Indices for 191228.8 dataset
indices= ...
[ 71 133
335 400
803 869
1606 1672
1739 1805
1943 2006
2475 2541
2610 2675 ];

nSets=size(indices,1);

if writeMicrographs

    CheckAndMakeDir(mergePath,1);
end;


for i=1:nSets
%     setPath=['Info_' sprintf('%04d',indices(i,1)) '/'];
    setPath='Info_redo/';
    CheckAndMakeDir(setPath,1);
    for k=indices(i,1):indices(i,2)
        [pa,nm,ex]=fileparts(allNames{k});
        mName_s=[nm(1:end-2) 'ms.mrc'];
        mName=[nm(1:end-2) 'm.mrc'];
        str1=['cp ' pa '/' nm '* ' setPath];
        str2=['cp ' mPath mName_s ' ' mergePath];
        str3=['cp ' mPath mName ' ' mergePath];
        if doWrite
            system(str1);
            if writeMicrographs
                system(str2);
                system(str3);
            end;
        end;
        disp(str1);
        if writeMicrographs
            disp(str2);
        end;
    end;
end;
return



%% Unsort
doWrite=1;
newInfoDir='Info/';
if exist(newInfoDir,'dir')
    error('Directory exists.');
end;
if doWrite
    CheckAndMakeDir(newInfoDir,1);
end;
d=dir;
for i=1:numel(d)
    if d(i).isdir && strncmp(d(i).name,'Info__',6)
        str=['cp ' d(i).name '/* ' newInfoDir];
        disp(str);
        if doWrite
            system(str);
        end;
    end;
end;