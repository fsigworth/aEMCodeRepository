function s=SetOptionValues(s,options,warningsOn)
% if options has a field with the same name as a field of s, overwrite
% s.field with options.field.
% Use this with default values in s, new values in options.
% By default, warnings are give in a field in options is not present in s.

if nargin<3
    warningsOn=1;
end;

sFields=fieldnames(s);
optFields=fieldnames(options);
for i=1:numel(optFields)
    p=strcmp(optFields{i},sFields);
    if any(p)
        ptr=find(p,1);
        s.(sFields{ptr})=options.(optFields{i});
    elseif warningsOn
        disp(['  Unrecognized option: ' optFields{i}]);
    end;
end;
