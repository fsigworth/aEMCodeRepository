function vals=miExtractFieldValues(allMis,fieldName,forceCells);
% From an allMis structure (cell array of structs), extract the given field
% of each struct into an array. A typical use is
% defoci=miExtractParameters(allMis,'ctf.defocus');
% Examples:
% fieldName = 'ctf' returns an nim x 1 cell array of structs.
% fieldName = 'vesicle.ok' typically returns a cell array (because the
%   numeric arrays are of variable size). Set the optional argument forceCells to
%   true to enforce a cell array.
% fieldName = 'ctf.defocus' returns a numeric nim x 1 vector.

maxDepth=3; % maximum levels of fields. That is, we can handle up to
% allMis{n}.a.b.c

if nargin<3
    forceCells=false;
end;

nim=numel(allMis);
fieldNames=strsplit(fieldName,'.');
vals=cell(nim,1);

for i=1:nim
    mi=allMis{i};
    nf=numel(fieldNames);
    ok=true;
    if isfield(mi,fieldNames{1})
        x=mi.(fieldNames{1});
        j=1;
        if nf>j && isfield(x,fieldNames{2})
            x=mi.(fieldNames{1}).(fieldNames{2});
            j=2;
            if nf>j && isfield(mi.(fieldNames{1}).(fieldNames{2}),fieldNames{3})
                x=mi.(fieldNames{1}).(fieldNames{2}).(fieldNames{3});
                if nf>maxDepth
                    j=maxDepth+1;
                    disp('Can''t handle fieldname depth >3.');
                    ok=false;
                    break
                end;
            end;
        end;
        vals{i}=x; % at least one field assigned.
    end;
end;
if ~ok
    vals=[];
    return
end;
if ~forceCells
    % See if we can return ordinary numeric values
    first=true;
    nvals=[];
    for i=1:nim
        if numel(vals{i})>0 && isnumeric(vals{i})
            if first
                sz=size(vals{i});
                nvals=zeros([nim sz]);
                first=false;
            else
                if ~all(size(vals{i})==sz)
                    break
                end;
            end;
            % convert to nmeric
            nvals(i,:)=vals{i}(:);
        end;
    end;
    if i==nim && numel(nvals)>0 % successful conversion
        vals=reshape(nvals,[nim sz]);
    end;
    % otherwise we leave vals as a cell array.
end;