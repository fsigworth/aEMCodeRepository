% CopyMiVesiclesAndParticles2.m
% Reads info from an old xxxmi.mat file into a new xxxmi.txt file


% Select the new info files from the file selector

disp('Selecting the new *mi.txt files');
[fname, pathName]=uigetfile('*mi.txt','Select new info files','multiselect','on');
if isnumeric(pathName) % File selection was cancelled
    return
end;
if ~iscell(fname)
    fname={fname};
end;
[rootPath, dirInfo]=ParsePath(pathName);
cd(pathName);  % Point to the mi file directory
%%
disp('Selecting the old *mi.mat files');
[fnameOld, pathOld]=uigetfile('*mi.mat','Select old info files','multiselect','on');
if isnumeric(pathOld) % File selection was cancelled
    return
end;
if ~iscell(fnameOld)
    fnameOld={fnameOld};
end;

%%
nOld=numel(fnameOld);
namesOld=cell(nOld,1);
for i=1:nOld
    nmOld=fnameOld(i);
    [pa,nm,ext]=fileparts(fnameOld{i});
    namesOld{i}=nm;
end;

nNew=numel(fname);
ptrsToOldNames=zeros(nNew,1);
for ind=1:nNew
    name=fname{ind};
    matches=strcmp(nm,namesOld);
    [pa,nm,ex]=fileparts(name);
    if any(matches)
        ptrsToOldNames(ind)=find(matches,1);
        disp([pathName name '  ' pathOld namesOld{ptrsToOldNames(ind)}]);
    else
        disp([pathName name '  --no match']);
        
    end;
end;
return

% 
% pAll=strfind(datesOld,dateStr);
% p=0;
% for j=1:nOld
%     if numel(pAll{j})>0
%         p=j;
%         break
%     end;
% end;
% if p>0
%     miNew=load(name);
%     mi=miNew.mi;
%     miOld=load([AddSlash(pathOld) fnameOld{p}]);
%     miOld=miOld.mi;
%     disp([name ' <--  ' fnameOld{p}]);
%     if numel(miOld.vesicle.x)>0
%         mi.mask=miOld.mask;
%         mi.vesicleModel=miOld.vesicleModel;
%         mi.vesicle=miOld.vesicle;
%         mi.boxSize=miOld.boxSize;
%         mi.particle=miOld.particle;
%         save(name,'mi');
%         return
%     end;
% end;
% 
% %%%
