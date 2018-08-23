function [name,ok]=CheckForImageOrZTiff(name)
% Check to see if an image with the given name exists.  If not, construct
% the compressed Tiff version of the name, and see if that exists.  If it
% does, return that name.  If not, ok is returned as 0, and the returned
% name is the same as the original one.
ok=1;
if ~exist(name,'file')
    %     Check to see if there is a compressed version
    [pa,nm,ex]=fileparts(name);
    altName{1}=[AddSlash(pa) nm 'z.tif'];  % compressed name
    altName{2}=[AddSlash(pa) nm '.z.tif'];  % other version
    altName{3}=[AddSlash(pa) nm '.tif'];   % was already tiff
    ok=0;
    for i=1:3
        if exist(altName{i},'file')
            name=altName{i};
            ok=1;
            break;
        end;
    end;
    if ~ok && nargout<2
        error(['File doesn''t exist: ' name]);
    end;
end;