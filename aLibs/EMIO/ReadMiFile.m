function [mi,actualName]=ReadMiFile(name,loadTextFirst,endFieldName)
% function [mi,actualName]=ReadMiFile(name,endFieldName)
% function [mi,actualName]=ReadMiFile(name,loadTextFirst,endFieldName)
% Default is to try to load the text version first.  loadTextFirst=false
% causes the .mat file version to be tried first.
% endFieldName is the name of the last field to be read in a truncated
% read (saves time in reading mi.txt files).
mi=[];
actualName='';
if nargin<1
    [name,pa]=uigetfile('*.txt');
    if isnumeric(pa)
        return;
    end;
    name=[pa name];
end;
    
if nargin<2
    loadTextFirst=true;
end;
if nargin<3
    if ischar(loadTextFirst)
        endFieldName=loadTextFirst;
        loadTextFirst=true;
    else
        endFieldName=char(26); % default is ctrl-Z
    end;
end;

[pa,nm,ex]=fileparts(name);
baseName=[AddSlash(pa) nm];
txtName=[baseName '.txt'];
matName=[baseName '.mat'];
mi=struct;
txtLoaded=false;
matLoaded=false;
if loadTextFirst && exist(txtName,'file')
    mi=ReadMiText(txtName,endFieldName);
    mi=ReadMiExtension(mi,txtName);  % partial read won't contain the mieFields, so no extension will be read.
    actualName=txtName;
    txtLoaded=true;
elseif exist(matName,'file')
    s=load(matName);
    mi=s.mi;
    actualName=matName;
    matLoaded=true;
else  % final attempt at text file
    if exist(txtName,'file')
        mi=ReadMiText(txtName,endFieldName);
        mi=ReadMiExtension(mi,txtName);
        actualName=txtName;
        txtLoaded=true;
    end;
end;
if ~(txtLoaded || matLoaded)
    error(['mi file could not be found: ' baseName]);
end

return

%% code to fix errors in writing, where spurious fields wound up in the
%  particle and vesicle sub-structures
nErrors=0;
mi0=meCreateMicrographInfoStruct14;
fields=fieldnames(mi0);  % get all the top-level fieldnames
structs={'vesicle' 'particle' 'vesicle'};
for j=1:numel(structs)
    sfield=structs{j};
    for i=1:numel(fields)
        if isfield(mi,sfield) && isfield(mi.(sfield),(fields{i}))
            fname=fields{i};
            doCopy=1;
            if isfield(mi,fname)  % there's also a root-level field.
                doCopy=numel(mi.(fname))<numel(mi.(sfield).(fname));
            end;
            if doCopy
                disp(['  ReadMiFile: Shifting the field: ' fname ' from ' sfield]);
                mi.(fname)=mi.(sfield).(fname);
            else
                disp(['  ReadMiFile: Removing the field: ' fname ' from ' sfield]);
            end;
            mi.(sfield)=rmfield(mi.(sfield),fname);
            nErrors=nErrors+1;
        end;
    end;
end;
if nErrors>0
    WriteMiFile(mi,actualName);
    disp(['  -corrected file written: ' actualName]);
end;
