function s=SetDefaultValues(defaults,s,checkFieldnames)
% function s=SetDefaultValues(defaults,s,checkFieldnames)
% Set default parameter values.  defaults and s are both structs.
% For each field that defaults has but s lacks, the default value is copied to s.
% if checkFieldnames is true (default is false) we list the fields of s
% that defaults doesn't have.

if nargin>2 && checkFieldnames  % we want warnings about unrecognized fields in s
    sNames=fieldnames(s);
    nfs=numel(sNames);
    for i=1:nfs
        if ~isfield(defaults,name)
            disp(['Unrecognized field: ' name]);
        end;
    end;
end;

% now copy over every value from defaults that s lacks.

names=fieldnames(defaults);
nf=numel(names);
for i=1:nf
    name=names{i};
    if ~isfield(s,name)
        s.(name)=defaults.(name);
    end;
end;

