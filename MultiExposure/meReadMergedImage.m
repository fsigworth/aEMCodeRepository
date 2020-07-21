function [m, mergeFullPath,ok]=meReadMergedImage(mi,doNorm,suffix)
% function [m mergeFullPath]=ReadMergedImage(mi,doNorm,suffix)
% Typical suffix values: '', 'v', 's', 'vs'
% if suffix(end)=='s' we first try to read a small image from mi.procPath_sm.
%
% Try to load the merged image m, and normalize it to total dose (if doNorm = 1;
% but note this is old-style normalization.  Now (2016) merged images are already
% normalized, so it's best to use the default which is no normalization).
% Compressed files ('z.tif') are no longer supported.
% If the image file isn't found, put up a file selector, and return
% both the image and the path where the image is found.
% m=0 if nothing is found.

if nargin<2
    doNorm=0;
end;
if nargin<3
    suffix='';
end;
mergeFullPath='';
m=0;
ok=false;
if isfield(mi,'procPath_sm')
    ppath_sm=mi.procPath_sm;
else
    ppath_sm=mi.procPath;
end;
for ibase=1:2 % Try mi.basePath, then pwd   
    ourSuffix=suffix;
    if ourSuffix(end)=='s'
        iname=[mi.basePath ppath_sm mi.baseFilename 'm' ourSuffix '.mrc'];
        if exist(iname,'file')
        ok=true;
        break;
        end;
        ourSuffix(end)='';
    end;
    if ~ok % look for a full-size file
        iname=[mi.basePath mi.procPath mi.baseFilename 'm' ourSuffix '.mrc'];
        if exist(iname,'file')
            ok=true;
            break;
        end;
    end;
    % for second round, try this path
    mi.basePath=AddSlash(pwd);
end; % for
if ok
    m=ReadEMFile(iname);
    mergeFullPath=iname;
elseif nargout<3  % if no error return given
    %%
    disp('Expected merged-image file not found:');
    disp(iname);
    [iname, mergeFullPath]=uigetfile('*m*.mrc','Find the merged image');
    if numel(iname)>1
        m=ReadEMFile([mergeFullPath iname]);
        mergeFullPath=[mergeFullPath '/' iname];
    end;
end;
if doNorm
    m=m/mi.doses(1);
end;
