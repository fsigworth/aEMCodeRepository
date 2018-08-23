function infoName=WriteMiFile(mi,name)
% function infoName=WriteMiFile(mi,writeText)
% function infoName=WriteMiFile(mi,name)
% assumes that we are already cd'd to the basePath.
% writeText is either a boolean, or a string containing
% the infoPath and filename to be used.
% The returned infoName is of the same form, InfoPath/Filename.txt,
% the same as used in the save operation.
if nargin<2
    name=1;
end;

if numel(name)>4 % actually a name
    infoName=name;
    [pa,nm,ext]=fileparts(infoName);
    q=strcmp(ext,{'.mat' '.txt'});
    if q(2)
        WriteMiText(mi,infoName);
    elseif q(1)
        save(infoName,'mi');
    else
        error(['Unrecognized extension on filename: ' ext]);
    end;
else  % create our own name, write as text
    if name
        infoName=[mi.infoPath mi.baseFilename 'mi.txt'];
        WriteMiText(mi,infoName);
    else
        infoName=[mi.infoPath mi.baseFilename 'mi.mat'];
        save(infoName,'mi');
    end;
end;
