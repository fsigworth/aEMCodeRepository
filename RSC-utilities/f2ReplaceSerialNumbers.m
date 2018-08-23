function f2ReplaceSerialNumbers(miNames)
% Update serial numbers in filenames e.g. from 3 to 4 digits
if nargin<1
    miNames={};
end;
serNoFormat='%04d';
doWrite=1;
recycleState=recycle;
recycle('on');

updateOnlyMerged=0;


if updateOnlyMerged
    pa=uigetdir;
    cd(pa);
    d=dir;
    for i=1:numel(d)
        name=d(i).name;
        if strndcmp(name,'z.tif',5)
            newName=ReplaceSerNo(name);
            if ~strcmp(newName,name)
            moveStr=['!mv ' name ' ' newName];
            if doWrite
                eval(moveStr);
            end;
            disp(moveStr);
            end;
        end;
    end;
else


nNames=numel(miNames);
for i=1:nNames
    oldMiName=miNames{i};
    mi=ReadMiFile(oldMiName);
    newMiName=ReplaceSerNo(oldMiName);
    oldBaseName=mi.baseFilename;
    mi.baseFilename=ReplaceSerNo(mi.baseFilename);
    if exist(mi.procPath,'dir') % we might have merged files
        oldMergeName=[mi.procPath oldBaseName 'm.mrc'];
        altMergeName=[mi.procPath oldBaseName 'mz.tif'];
        if ~exist(oldMergeName,'file')
            oldMergeName=altMergeName;
        end;
        if exist(oldMergeName,'file')
            newMergeName=ReplaceSerNo(oldMergeName);
            moveStr=['!mv ' oldMergeName ' ' newMergeName];
            if doWrite
                eval(moveStr);
            end;
            disp(['        ' moveStr]);
        end;
    end;
    %         Convert the micrograph names, if present
    if isfield(mi,'imageFilenames') && numel(mi.imageFilenames)>0
        imgFilenames=mi.imageFilenames;
        mi.imageFilenames={};
        for j=1:numel(imgFilenames)
            oldImageName=imgFilenames{j};
            [pao, nmo, exo]=fileparts(oldImageName);
            if ~exist([mi.imagePath oldImageName],'file')
                altOldImageName=[nmo 'z.tif'];
                if ~exist([mi.imagePath altOldImageName],'file')
                    error(['Can''t find the file ' mi.imagePath oldImageName]);
                else
                    oldImageName=altOldImageName;
                end;
            end;
            newImageName=ReplaceSerNo(oldImageName);
            mi.imageFilenames{j,1}=newImageName;
            
            moveStr=['!mv ' mi.imagePath oldImageName ' ' mi.imagePath newImageName];
            if doWrite
                eval(moveStr);
            end;
            disp(['        ' moveStr]);
        end;
    end;
    disp(['Writing: ' newMiName])
    if doWrite
        WriteMiFile(mi,newMiName);
        delete(oldMiName);
        WriteMiFile(mi,newMiName); % write again just to make sure.
    end;
end;
end;
recycle(recycleState);

    function newName=ReplaceSerNo(oldName)
        [pa,nm,ex]=fileparts(oldName);
        ptrs=strfind(nm,'_');
        if numel(ptrs)<2
            error(['bad filename: ' nm]);
        end;
        serNo=str2double(nm(ptrs(1)+1:ptrs(2)-1));
        newNm=[nm(1:ptrs(1)) sprintf(serNoFormat,serNo) nm(ptrs(2):end)];
        newName=[AddSlash(pa) newNm ex];
    end
end
