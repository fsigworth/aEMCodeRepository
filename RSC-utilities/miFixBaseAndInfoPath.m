% miFixBaseAndInfoPath
% The user selects some mi files.  The base directory name (the one above where
% the files are now) is written into the mi.basePath field, and the infoPath
% field is likewise updated.

% Have the user select some mi files
[fname pa]=uigetfile('*mi.mat','Select mi files to fix','multiselect','on');
[rootPath infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;
cd(rootPath);


nfiles=numel(fname)
if nfiles>0
    disp('Changing the mi.basePath to:');
    disp(['    ' rootPath]);
    disp('and mi.infoPath to:');
    disp(['    ' infoPath]);
    ok=input('ok [y]? ','s');
    if numel(ok)<1 || lower(ok(1))=='y'  % default is 'yes'
        
        for fileIndex=1:nfiles; % Operate on a single mi file
            miName=[infoPath fname{fileIndex}];
            disp(['Changing ' miName]);
            load(miName);
            mi.basePath=rootPath;
            mi.infoPath=infoPath;
            save(miName,'mi');
        end;
    end;
else
    disp('No files selected.');
end;
disp('Done.');
