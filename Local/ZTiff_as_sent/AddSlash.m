function pathstr=AddSlash(pathstr)
% function pathstr=AddSlash(pathstr)
% given a path string, check that it ends with a the 
% platform's fileseperator (e.g. slash).  If not, add one.  This allows the
% path strings returned by other functions to be used in constructing full
% filenames, e.g.
% fullName=[AddSlash(pathName) filename];

n=numel(pathstr);
if n>0 && pathstr(n) ~=filesep
    pathstr=[pathstr filesep];
end;
