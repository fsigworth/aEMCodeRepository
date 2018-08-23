function s=SetDefaultValues(s,options)
% if options has a field with the same name as a field of s, overwrite
% s.field with options.field.
% Use this with default values set to s, new values to options.
    sFields=getfields(s);
    optFields=getFields(options);
    for i=1:numel(sFields)
        p=strcmp(sFields{i},optFields);
        if any(p)
            s.(sFields(i))=options.(sFields{i});
        end;
    end;