function s=AddStructFields(s,s1)
% Add the numeric value of each field of s1 to the corresponding field of
% s, and return the result.

names=fieldnames(s1);
for i=1:numel(names)
    nm=names{i};
    if ~isfield(s,nm)
        error(['Structure s lacks the field ' nm]);
    end;
    s.(nm)=s.(nm)+s1.(nm);
end;


