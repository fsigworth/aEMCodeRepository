function [is,js]=reGetFlagEntries(emFlags)
switch emFlags.type
    case 'logical'
        sz=size(emFlags.flags);
        nRows=prod(sz(1:3));  % nalphas x nrefs x nvols
        nCols=sz(4);  % nImgs
        indices=find(emFlags.flags);
        [is,js]=ind2sub([nRows nCols],indices);
        is=uint32(is);
        js=uint32(js);
    otherwise
        disp(['Unknown flag type: ' emFlags.type]);
end;

end
