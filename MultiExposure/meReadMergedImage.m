function [m, mergeFullPath,ok]=meReadMergedImage(mi,doNorm,suffix)
% function [m mergeFullPath]=ReadMergedImage(mi,doNorm,suffix)
% Try to load the merged image m, and normalize it to total dose (if doNorm = 1;
% but note this is old-style normalization.  Now (2016 merged images are already
% normalized, so it's best to use the default which is no normalization).
% The usual filename is of the form *m.mrc, but if suffix='v' it is *mv.mrc
% Alternatively, it may be *mz.tif or *mvz.tif in these cases, for
% compressed files.
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

iname=[mi.basePath mi.procPath mi.baseFilename 'm' suffix '.mrc'];
[iname,ok]=CheckForImageOrZTiff(iname);
if ~ok % try for the present path instead
    mi.basePath=AddSlash(pwd);
    iname=[mi.basePath mi.procPath mi.baseFilename 'm' suffix '.mrc'];
    [iname,ok]=CheckForImageOrZTiff(iname);
end;

if ok
    m=ReadEMFile(iname);
    mergeFullPath=iname;
elseif nargout<3  % if no error return given
    %%
    disp('Expected merged-image file not found:');
    disp(iname);
    [iname, mergeFullPath]=uigetfile(['*m' suffix '*'],'Find the merged image');
    if numel(iname)>1
        m=ReadEMFile([mergeFullPath iname]);
        mergeFullPath=[mergeFullPath '/' iname];
    end;
end;
if doNorm
    m=m/mi.doses(1);
end;
