% miChangeField
% The user selects one or more mi files.  The user selects a field to
% change, and gives the new value.  It is replaced in each mi file.

% Have the user select some mi files
[fname, pa]=uigetfile('*mi.mat','Select mi files to change','multiselect','on');
[rootPath, infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;
cd(rootPath);

nFiles=numel(fname);
if nFiles<1
    return
elseif nFiles==1
    fileStr=' file';
else
    fileStr=' files';
end;

disp([num2str(nFiles) fileStr ' starting with ' fname{1}]);
ok=1;
while ok
    miName=[infoPath fname{1}];
    load(miName);
    disp(mi);
    field=input('Name of field to change, or q to quit ? ','s');
    if lower(field(1))=='q'
        disp('Done.');
        return
    end;
    if isfield(mi,field)  % is a valid field
        valOk=1;
        val=mi.(field);
        if isa(val,'numeric') && size(val,1)<2
            val=MyInput('new value',val);
        elseif isa(val,'char') && size(val,1)<2
            oldVal=val;
            if strcmp(field,'basePath')  % special case
                   q=input('Reset to current directory [y]?','s');
                   if numel(q)<1 || lower(q)=='y'
                     val=rootPath;
                   else
                     val=input(['new value [' val ']? '],'s');  
                   end;
            else
                val=input(['new value [' val ']? '],'s');
            end;
            if numel(val)<1
                val=oldVal;
            end;
        else
            disp(' Can''t enter this field, sorry.');
            valOk=0;
        end;
        if valOk
            q=input('Are you sure [y]?','s');
            if numel(q)<1 || lower(q)=='y'
                disp(['Changing ' num2str(nFiles) fileStr]);
                for fileIndex=1:nFiles; % Operate on a single mi file
                    miName=[infoPath fname{fileIndex}];
                    load(miName);
                    mi.(field)=val;
                    save(miName,'mi');
                    disp(miName);
                end;
            end;
        end;
    else
        disp('Not a valid field name.');
    end;
end;
disp('Done.');
