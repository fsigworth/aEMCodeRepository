function [names,pathName,ok]=uiGetFilenames(filter,windowTitle)
% Friendly version of uigetfile.  Always returns names as a cell array;
% path has a slash at the end.  ok is false and the other returned values
% are empty if cancel is clicked. The function can be called with no
% arguments.
% Filter is a cell array of filter strings, e.g. {'*mi.mat' '*mi.txt'};
% default is no filter.
if nargin<1
    filter='*.*';
end;
if nargin<2
    windowTitle='';
end;
[names, pathName,ok]=uigetfile(filter,windowTitle,'multiselect','on');
if isnumeric(pathName)
    ok=false;
    names=cell(0);
    pathName='';
    return
end;
if ~iscell(names)
    names={names};
end;
ok=true;
return