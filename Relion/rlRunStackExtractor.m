% rlRunStackExtractor

load workDirs4.mat
siNames=cell(numel(workDirs),1);
for idir=1:numel(workDirs)
    disp(['---working on ' workDirs{idir} '---']);
    cd(workDirs{idir});
    % StackExtractor3;
    
    siPath='Stack2/';
    d=dir(siPath);
    siName='';
    for j=1:numel(d)
        if strndcmp(d(j).name,'tsi.mat')
            siName=d(j).name;
            break;
        end;
    end;
    disp(siName);
    siNames{idir}=siName;
    if numel(siName)>5
        SiToStarFile;
    else
        disp(['No tsi.mat file found in ' AddSlash(pwd) siPath]);
    end;
end;

return
%% Merge the star files


