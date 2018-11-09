function s=SetDefaultValues(defaults,s)
% function s=SetDefaultValues(defaults,s)
% Set default parameter values.  defaults and s are both structs.
% For each field that defaults has but s lacks, the default value is copied to s.

names=fieldnames(defaults);    
    nf=numel(names);
    for i=1:nf
        name=names{i};
        if ~isfield(s,name)
            s.(name)=defaults.(name);
        end;
    end;
