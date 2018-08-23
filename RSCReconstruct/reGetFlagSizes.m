function [nAlphas, nRefs, nVols, nImgs, is, js]=reGetFlagSizes(emFlags)
switch emFlags.type
    case 'logical'
        sz=size(emFlags.flags);
        nAlphas=sz(1);
        nRefs=sz(2);
        nVols=sz(3);
        nImgs=sz(4);
        
        if nargout>4
            nRows=prod(sz(1:3));  % nalphas x nrefs x nvols
            nCols=sz(4);  % nImgs
            indices=find(emFlags.flags);
            [is,js]=ind2sub([nRows nCols],indices);
            is=uint32(is);
            js=uint32(js);
        end;
    otherwise
        disp(['Unknown flag type: ' emFlags.type]);
end;
end
