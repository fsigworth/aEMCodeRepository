% MiLoadAll.m
% Use the file selector to choose a set of mi files, or look in the 
% current directory.  Load all of them into
% a single cell array, and store all their names too.

% Use the later sections of code to find defocus values, and sort micrographs
% according to the defocus values.

getNamesOnly=0;
batchMode=1;
if batchMode
    outPath='';
    infoPath='Info/';
    names=f2FindInfoFiles(infoPath);
else
    [names, pathName]=uigetfile({'*mi.txt';'*mi.mat'},'Select mi files','multiselect','on');
    if isnumeric(pathName) % File selection was cancelled
        return
    end;
    if ~iscell(names)
        names={names};
    end;
    [rootPath, infoPath]=ParsePath(pathName);
    cd(rootPath);
    outPath=infoPath;
end;
%%
allMis=cell(0);
nNames=numel(names);
disp([num2str(nNames) ' mi files.']);
allNames=cell(nNames,1);
for i=1:nNames
    disp(names{i});
    %     allNames{i}=[infoPath names{i}];
    allNames{i}=names{i};
    if ~getNamesOnly
        allMis{i,1}=ReadMiFile(allNames{i});
    end;
end;
save([outPath 'allNames.mat'],'allNames');
disp([outPath 'allNames.mat saved.']);
if ~getNamesOnly
    save([outPath 'allMis.mat'],'allMis','allNames');
    disp([outPath 'allMis.mat saved.']);
end;

return

%% Make histogram of defocus values
defs=zeros(nNames,1);
for i=1:nNames
    defs(i)=allMis{i}.ctf(1).defocus;
end;

hist(defs,50);

return

%% Find the files with defocus closest to the target value. Put the micrographs into
%% an array of images, or else copy the mrc files from a remote place.
defTarget=2.26;
numToLoad=20;
jpegDir='jpeg/';
fc=.2;
nsds=4;
doWrite=0;
% CheckAndMakeDir(jpegDir,1);

copyRemoteMrcs=1;

if copyRemoteMrcs
    remotePath='/Volumes/Drobo4/200922/Merged/';
    % localPath='/Volumes/D257/Hideki/200922/Merged/';
    localPath='Merged/';
    h=fopen('CopyScript.sh','w');
end

[sDeltas,sInds]=sort(abs(defTarget-defs));
foundMis=allMis(sInds);
foundDefs=defs(sInds);
foundNames=allNames(sInds);

suffix={'s' 'vs'};
suffixFull={'m' 'mv'};

for i=1:numToLoad
    mi=foundMis{i};
    for j=1:2
        if copyRemoteMrcs
            
            fileName=[mi.baseFilename suffixFull{j} '.mrc'];
            str=['cp ' remotePath fileName ' ' localPath fileName];
             disp(str);
            fprintf(h,'%s\n',str);
            
            
        else  % read the 'ms' files locally, create the ms image array
%             and make jpegs
            [m, mergeFullName,ok]=meReadMergedImage(mi,0,suffix{j});
            disp(mergeFullName);
            if ~ok
                return
            end;
            if i==1 && j==1
                ms=zeros([size(m) numToLoad 2],'single');
            end;
            ms(:,:,i,j)=m;
            mscl=uint8(imscale(GaussFilt(m,fc),256,nsds));
            imaga(mscl);
            jName=[jpegDir num2str(i) '_' mi.baseFilename 'm' suffix{j} '.jpg'];
            title([num2str(mi.ctf(1).defocus,3) '  ' jName], 'interpreter', 'none');
            drawnow;
            if doWrite
                imwrite(mscl,jName);
            end;
            
        end;
    end
end;
if copyRemoteMrcs
    fclose(h);
    save DefSorted.mat foundDefs foundMis foundNames
else
    save DefSorted.mat foundDefs foundMis foundNames ms
end;
