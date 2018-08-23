% RenameFiles.m
monthStr='Jan';

        % Select info files from the file selector
        [fnameOld, pathOld]=uigetfile('*rscc.mat','Select old info files','multiselect','on');
        if isnumeric(pathOld) % File selection was cancelled
            return
        end;
        if ~iscell(fnameOld)
            fname={fnameOld};
        end;
        cd(pathOld);
        
% Get all the old file dates
nOld=numel(fnameOld);

for ind=1:nOld
    name=fnameOld{ind};
    p=strfind(name,monthStr);
    baseName=name(1:p+13);
    newName=[baseName 'rscc.mat'];
    rscc=load(name);
     save(newName,'-struct','rscc');
    disp([name '    ' newName]);
end;
