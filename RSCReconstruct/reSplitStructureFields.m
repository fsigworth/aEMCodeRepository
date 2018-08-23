function s1=reSplitStructureFields(s,flags)
% Sample the structure s to yield the smaller structure s1.  This is used
% in dividing up the reconstruction problem among groups and slices.
% Extensive fields are sampled according to the logical vector flags.
% Intensive fields are simply copied into the output structure.  The names
% of the intensive fields are found in the cell array s.intensiveFields;
% also fields having a final dimension < numel(flags) is assumed to be
% intensive.

warningOn=0;  % If set, warn when an extensive field has a final dimension 
              % that is too small and is being treated as intensive.

if isfield(s,'intensiveFields')
    intFields=s.intensiveFields;
    intFields(end+1)={'intensiveFields'};  % add this one too.
else
    intFields={};
end;
allFields=fieldnames(s);  % default extensive field list
extFields=allFields;

% Remove intensive fields from the list.
intensOk=false(numel(intFields),1);
for i=1:numel(intFields)
    q=strcmp(extFields,intFields{i});
    if any(q)
        extFields(q)=[];
        intensOk(i)=true;
    end;
end;
% Prune any unused intensive fields
intFields=intFields(intensOk);

s1=struct;
flags=logical(flags);
for i=1:numel(extFields)
    name=extFields{i};
    sz=size(s.(name));
    if sz(end)==1  % remove trailing singleton dimension
        sz(end)=[];
    end;
    ndims=numel(sz);
    if sz(end)<numel(flags) % too small to be extensive
        if warningOn
            warning(['the field ' name ' assumed to be intensive']);
        end;
        intFields{end+1}=name;
    else
        switch ndims
            case 1
                s1.(name)=s.(name)(flags,1);
            case 2
                s1.(name)=s.(name)(:,flags);
            case 3
                s1.(name)=s.(name)(:,:,flags);
            case 4
                s1.(name)=s.(name)(:,:,:,flags);
            otherwise
                error(['Too many dimensions: ' num2str(ndims)]);
        end;
    end;
end;

for i=1:numel(intFields)
    name=intFields{i};
    s1.(name)=s.(name);
end;
