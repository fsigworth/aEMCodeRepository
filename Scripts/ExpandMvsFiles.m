% PadSmallImages
% Fix the buggy outout of rsRefineVesicleFits when the mvs.mrc file is the
% wrong size.
load allMis;
nim=numel(allMis);

for i=1:nim
    mi=allMis{i}; 
    mvsName=[mi.procPath_sm mi.baseFilename 'mvs.mrc'];
    if exist(mvsName,'file')
        [mvs,s]=ReadMRC(mvsName);
        disp(mvsName);
        WriteMRC(Crop(mvs,960),s.pixA,[mvsName 'x']);
    end;
end;

return
%%
doExec=1;
for i=1:nim
    mi=allMis{i}; 
    mvsName=[mi.procPath_sm mi.baseFilename 'mvs.mrc'];
    mvsxName=[mvsName 'x'];
    if exist(mvsxName,'file')
        str=['mv ' mvsxName ' ' mvsName];
        if doExec
            system(str);
        end;
        disp(str);
    end;
end;
