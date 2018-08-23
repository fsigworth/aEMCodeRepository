% rlSetActiveFlagsFromStar

disp('Getting a si.mat file from the original images');
[siName, siPath]=uigetfile('*si.mat','Select stack info file');
if isnumeric(siPath)  % user has clicked Cancel
    return
end;

load([AddSlash(siPath) siName]);
ShowActives(si);

goOn=1;

while goOn
    disp('Getting a data.star file from the classification.');
    [stName, stPath]=uigetfile('*.star','Select data star file');
    if isnumeric(stPath)  % user has clicked Cancel
        return
    else
        cd(stPath);
    end;
    %%
    disp(['Reading ' stName '...']);
    [stBlocks,stData,ok]=ReadStarFile(stName);
    disp('done.');
    %%
    dat=stData{1};
    nim=numel(dat.rlnDefocusU);
    if isfield(dat,'rlnClassNumber')
        cls=dat.rlnClassNumber;
    else
        cls=zeros(nim,1);
    end;
    
    ind=zeros(nim,1); % original image index
    active=false(size(si.activeFlags,1),1);
    for i=1:nim
        ind(i)=rlDecodeImageName(dat.rlnImageName{i});
        active(ind(i))=true;
    end;
    disp([num2str(sum(active)) ' active images']);
    txt=input('Flag log text (blank to skip): ','s');
    if numel(txt)>1
        si.activeFlags(:,end+1)=active;
        si.activeFlagLog{end+1,1}=txt;
    end;
    q=input('Add more? ','s');
    goOn=(lower(q)=='y');
end;

ShowActives(si);

q=input('Overwrite the si file? ','s');
if lower(q)=='y'
    outSiPath=siPath;
    outSiName=siName;
else
    
    [outSiName, outSiPath]=uiputfile('*si.mat','Stack info file to write');
    if isnumeric(outSiPath)
        return
    end;
end;
disp(['Writing ' AddSlash(outSiPath) outSiName]);
save([AddSlash(outSiPath) outSiName]);
disp('done.');


function ShowActives(si)
for i=1:size(si.activeFlags,2)
    disp([num2str(i) '  ' si.activeFlagLog{i} '  ' num2str(sum(si.activeFlags(:,i)))]);
end;
end