function pathstr=AddSlash(pathstr)
% function pathstr=AddSlash(pathstr)
% given a path string, check that it ends with a the 
% platform's fileseperator (e.g. slash).  If not, add one.  This allows the
% path strings returned by other functions to be used in constructing fill
% filenames, e.g.
% fullName=[AddSlash(pathName) filename];
% This function also searches for and replaces improper file separators.
% Specifically it replaces all forward slashes with backslashes in the
% Windows environment, and vice versa for Unix environments.

% Simple implementation, assumes 1 byte = 1 character
% pathStrings='/\';
% if numel(pathstr)>0 && ~any(pathstr(end) == pathStrings)
%     pathstr=[pathstr filesep];
% end;
% 
% badSep=pathStrings(pathStrings~=filesep);
% q=strfind(pathStrings,badSep);
% for i=q
%     pathstr(i)=filesep;
% end;
% return

% Fancier implementation using cells and string operations.
pathStrings={'/' '\'};
if numel(pathstr)>0  % strcmp doesn't like empty strings
    endSeps=strcmp(pathstr(end),pathStrings);
    if ~any(endSeps)
        pathstr=[pathstr filesep];
    end;
end;

badSep=pathStrings{~strcmp(pathStrings,filesep)};
q=strfind(pathstr,badSep);  % This will work if there are only two possibilities.
for i=q
    pathstr(i)=filesep;
end;
