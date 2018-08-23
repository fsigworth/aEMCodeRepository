function [m name]=meGetPosterImage(mi)
% function [m name]=meGetPosterImage(mi)
% Given the micrograph info structure mi, get a representative image to
% show.  We first try to find a merged image *m.jpg or *m.mrc.  If that
% fails, we load the largest-defocus raw image.  The filename of the file
% that was used is returned.  m=0 if nothing found.

if mi.version<10  % old version doesn't have paths
    disp(['mi has an old version: ' num2str(mi.version)]);
    m=0;
    return
end;

% First search for a merged image.
mergedNameJ=[mi.procPath mi.baseFilename 'm.jpg'];
mergedNameM=[mi.procPath mi.baseFilename 'm.mrc'];
rawName=[mi.imagePath mi.imageFilenames{numel(mi.imageFilenames)}];
if FileExists(mergedNameJ)
    name=mergedNameJ;
    m=single(imread(name));
elseif FileExists(mergedNameM)
    name=mergedNameM;
    m=ReadEMFile(name);
elseif FileExists(rawName)
    name=rawName;
    m=ReadEMFile(name); % get the last of the raw images
    me=mean(m(:));
    sd=std(m(:));
    m(abs(m-me)>6*sd)=me;  % remove outliers
else
    m=0;
    name='';
    disp(['Image files not found: ' mergedNameM ' | ' rawName]);
end;
