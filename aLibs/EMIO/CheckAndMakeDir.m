function alreadyExists=CheckAndMakeDir(path,doDisplay)
% function alreadyExists=CheckAndMakeDir(path,doDisplay)
% Returns alreadyExists=true if the path already exists; otherwise, makes a new
% directory.  If doDisplay is true, it also writes out a message if
% a directory has been created.  (Default is doDisplay=false.)
if nargin<2
    doDisplay=false;
end;
    
alreadyExists=numel(path)<1 || exist(path,'dir');
if ~alreadyExists
    mkdir(path);
    if doDisplay
        disp(['Created the directory ' path]);
    end;
end;
