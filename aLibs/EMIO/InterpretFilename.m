function name=InterpretFilename(input)
% Allows a filename to be passed as a directory structure (use the 'name'
% field) or as a cell array element.
if isa(input,'char')
    name=input;
    return;
elseif isa(input,'struct')
    if isfield(input,'name')
        name=input.name;
    else
        error('Unrecognized structure passed as filename.');
    end;
elseif isa(input,'cell')
    name=input{1};
    if ~isa(name,'char')
        error('Unrecognized cell array passed as filename.');
    end;
else
    error(['Unrecognized class for filename: ' class(input)]);
end;
