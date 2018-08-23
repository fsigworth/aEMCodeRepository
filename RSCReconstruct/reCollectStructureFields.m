function so=reCollectStructureFields(si,flags)
% function so=reCollectStructureFields(si,flags)
% Collect the elements from a cell array of structures si.  For extensive fields we look at
% si{1} to obtain fieldnames and number of dimensions.  We copy along the last
% non-singleton dimension, and the flags tell us which elements of so come
% each element of si.  For example, we copy the elements of the 2d array
% si{i}.(name) into so.(name)(:,flags(:,i)).
% Each intensive field is simply summed to yield the field in so.  The names
% of intensive fields is taken from the cell array si{1}.intensiveFields.

% Example:
% si=cell(2,1);
% m.a=[1;2;3];
% m.b=[0 1 2;1 2 3];
% m.c=1:10;
% m.intensiveFields={'c'};
% si={m m};
% flags=logical([1 0; 1 0; 1 0; 0 1; 0 1; 0 1]);
% %
% so=reCollectStructureFields(si,flags)
%   Result:
% so.a' =
%      1     2     3     1     2     3
% so.b =
%      0     1     2     0     1     2
%      1     2     3     1     2     3
% so.c = 120
warningsOn=0;

if isfield(si{1},'intensiveFields')
    intFields=si{1}.intensiveFields;
else
    intFields={};
end;
    intFieldsi=intFields;  % include the "intensiveFields" field.
    intFieldsi(end+1)={'intensiveFields'};

[n,nSlices]=size(flags);
allFields=fieldnames(si{1});
extFields=allFields;

for i=1:numel(intFieldsi)
    q=strcmp(extFields,intFieldsi{i});
    if any(q)
        extFields(q)=[];
    end;
end;

flags=logical(flags);  % has to be logical to allow indexing.
so=struct;
for i=1:numel(extFields)
    name=extFields{i};
    szi=size(si{1}.(name));
    ndims=numel(szi);
    if ndims>1 && szi(ndims)==1  % remove trailing singleton dimension
        ndims=ndims-1;
    end;
    if sum(flags(:,1),1)~=szi(ndims)
        if warningsOn
            warning(['Size of the extensive field ' name ' doesn''t match flags']);
        end;
        intFields{end+1}=name;
    else
        szo=szi;
        szo(ndims)=n;  % make the last output dimension equal to n
        if isa(si{1}.(name),'logical')
            t=false(szo);
        else        
            t=zeros(szo,'like',si{1}.(name));
        end;
        for j=1:nSlices
            switch ndims
                case 1
                    t(flags(:,j))=si{j}.(name);
                case 2
                    t(:,flags(:,j))=si{j}.(name);
                case 3
                    t(:,:,flags(:,j))=si{j}.(name);
                case 4
                    t(:,:,:,flags(:,j))=si{j}.(name);
                otherwise
                    error(['Too many dimensions: ' num2str(ndims)]);
            end;
        end;
        so.(name)=t;
    end;
end;
for i=1:numel(intFields)
    name=intFields{i};
    sz=size(si{1}.(name));
    so.(name)=si{1}.(name);
    for j=2:nSlices
        so.(name)=so.(name)+si{j}.(name);
    end;
end;
so.intensiveFields=si{1}.intensiveFields;
