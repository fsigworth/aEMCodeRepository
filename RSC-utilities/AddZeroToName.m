% AddZeroToName
% special util for 160909/KvLipo121_2 info files.
% modify mi file names to add a leading zero.  Modify the baseFilename too,
% and modify all the data files.
% Assume we're in the Info directory.

doIt=1;
% cd('/Volumes/louise scratch/transfer/160909/KvLipo121_2')
cd('/Volumes/D215/160909/KvLipo121_2')

infoDir='Info/';
dirNames={'Micrograph/'; 'Merged/'; 'Merged/jpeg/';...
    'Merged/filtered/'; 'Jpeg/'};

%   Pick up the various directory contents
nDirs=numel(dirNames);
dList=cell(nDirs,1);
for i=1:nDirs
    d=dir(dirNames{i});
    nms=cell(numel(d),1);
    for j=1:numel(d)
        nms{j}=d(j).name;
    end;
    dList{i}=nms;
end;


usPtr=2;  % underscore preceding the index string
oldNumSize=3;  % number of digits
newNumSize=4;

dInfo=dir(infoDir);
ndi=numel(dInfo);
iNames=cell(0,1);
iSizes=zeros(0,1);
k=0;
for i=1:ndi  % scan for suitable filenames: has underscores and ends with 'mi.txt'
    name=dInfo(i).name;
    u=strfind(name,'_');
    if numel(name)>6 && strcmp(name(end-5:end),'mi.txt')...
            && numel(u)>usPtr % it's an mi file, has underscores
        k=k+1;
        %         disp(name)
        iNames{k,1}=name;
        iSizes(k,1)=dInfo(i).bytes;
    end;
end;
ndi=k;

% disp(' ');

toDel=cell(0,1);
numDel=0;

for i=1:ndi
    name=iNames{i};
    u=strfind(name,'_');
    disp(' ');
    disp(name);
    if u(usPtr+1)-u(usPtr)-1==oldNumSize  % found a file of the old size
        %             See if there's another info file with a leading zero
        xName=[name(1:u(usPtr)) '0' name(u(usPtr)+1:end)];
        disp(['Extended name: ' xName]);
        found=strcmp(xName,iNames);
        if any(found)
            j=find(found,1);  % j is index of longer file.
            if iSizes(j) > iSizes(i) % lengthened file is larger
                numDel=numDel+1;
                toDel{numDel,1}=name;  % file to delete
                name=xName;  % mi file to read
            else
                numDel=numDel+1;
                toDel{numDel,1}=xName;  % delete the short one
            end;
            disp(['** two mi files ' name ' -> ' xName ' ' num2str([iSizes(i) iSizes(j)])])
        end;
        for id=1:numDel
            str=['! rm ' infoDir toDel{id}];
            disp(str);
            if doIt, eval(str); end;
        end;
        numDel=0;
    end;  % handle short filenames
    disp(['Reading ' name])
    %         Now check the baseFilename
    mi=ReadMiFile([infoDir name]);
    mi0=mi;
    %%
    u=strfind(mi0.baseFilename,'_');
    if u(usPtr+1)-u(usPtr)-1==oldNumSize  % needs expansion
        disp(['mi filename is ' name]);
        mi.baseFilename=[mi0.baseFilename(1:u(usPtr)) '0' mi0.baseFilename(u(usPtr)+1:end)];
        disp(['baseFilename: ' mi0.baseFilename ' >> ' mi.baseFilename]);
        %         Fix the imageFilenames
        for im=1:numel(mi.imageFilenames)
            nm=mi.imageFilenames{im};
            u=strfind(nm,'_');
            if u(usPtr+1)-u(usPtr)-1==oldNumSize  % needs expansion
                mi.imageFilenames{im}=[nm(1:u(usPtr)) '0' nm(u(usPtr)+1:end)];
            end;
        end;
        disp(mi.imageFilenames);
        %             Done modifyin the mi file
        str=['!rm ' infoDir name];  % remove the old Mi file
        disp(str);
        disp('Write mi file');
        if doIt, eval(str); end;
        if doIt, WriteMiFile(mi,1); end;
        
        %        Now rename all the other files
        baseName=mi0.baseFilename;
        baseLen=numel(baseName);
        for j=1:nDirs
            q=find(strncmp(baseName,dList{j},baseLen));
            for k=1:numel(q)
                oldName=dList{j}{q(k)};
                oldPart=oldName(baseLen+1:end);
                fullOldName=[dirNames{j} oldName];
                fullNewName=[dirNames{j} mi.baseFilename oldPart];
                if ~strcmp(fullOldName, fullNewName)
                    str=['!mv ' fullOldName ' ' fullNewName];
                    disp(str);
                    if doIt, eval(str); end;
                end;
            end;
        end;
        
        
    end;
end;
%
%
% ex1=d(10).name;
% nx1=numel(ex1);
%
% zeroPos=8;
%
% ex2=d(120).name;
% nx2=numel(ex2);
%
% nd=numel(d);
% changed=false(nd,1);
% for i=1:nd
%     nm=d(i).name;
%     if numel(nm)==nx1  % a short name
%         ind=str2double(nm(zeroPos:zeroPos+2));
%         if ind<=493
%             changed(i)=true;
%             nm2=[nm(1:zeroPos-1) '0' nm(zeroPos:end)];
%             str=['!mv ' nm ' ' nm2];
%             disp(str);
%             if doWrite
%                 eval(str);
%             end;
%         end;
%     end;
% end;
