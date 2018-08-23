function flags=rlMapActiveFlags(selection,selFlags,default)
% function flags=rlMapActiveFlags(selection,selFlags,default)
% Suppose out of a set of particles img(:,:,selection) active particles are
% found, marked with selFlags.  Map these back to the full set of particles
% img(:,:,:).  All arguments are boolean vectors.
% default (normally false) is the flags value(s) when unset by the
% selection.
if nargin<3
    default=false;
end;
flags=false(size(selection)) | default;
selInds=find(selection);
flags(selInds(selFlags))=true;
