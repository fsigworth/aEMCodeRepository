% miChangeField
% The user selects one or more mi files.  The user selects a field to
% change, and gives the new value.  It is replaced in each mi file.

% % Have the user select some mi files
[fname, pa]=uigetfile('*mi.txt','Select mi files to change','multiselect','on');
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
miName=[infoPath fname{1}];
mi=ReadMiFile(miName);
disp(mi);

%%
fields=cell(0,1);
vals=cell(0,1);
nf=0; % number of fields to change
selecting=1;
while selecting
    field=input('Name of field to change or q to proceed ? ','s');
    if lower(field(1))=='q'
        selecting=0;
    else
        if isfield(mi,field) % is a valid field
            valOk=1;
            val=mi.(field);
        else
            field=input('new field name ?', 's');
            if numel(field)>0
                valOk=1;
                val='   ';
            else
                valOk=0;
            end;
        end;
        if valOk
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
                nf=nf+1;
                fields{nf}=field;
                vals{nf}=val;
                disp([num2str(nf) ' fields.']);
            end;
        elseif lower(field(1))~='q'
            disp('Not a field, try again');
        end;
    end;
end; % while
fields
vals
if nf<1
    disp('No fields selected. Done.');
    return
end;
%%
q=input(['Ready to change ' num2str(nFiles) ' mi files [y]?'],'s');
if numel(q)<1 || lower(q)=='y'
    disp(['Changing ' num2str(nFiles) fileStr]);
    for fileIndex=1:nFiles; % Operate on a single mi file
        miName=[infoPath fname{fileIndex}];
        mi=ReadMiFile(miName);
        for i=1:nf
            mi.(fields{i})=vals{i};
        end;
        WriteMiFile(mi,miName);
        disp([num2str(fileIndex) ': ' miName]);
    end;
end;

disp('Done.');
