function mi=meInsertMask(msk,mi,index,mode)
% function mi=meInsertMask(msk,mi,index,mode)
% Given a binary image msk (masked points are zeros), encode it and insert
% it into the mi structure as mi.mask(index).
% This information is decoded by meGetMask.

if nargin<4
    mode='AND';
end;
if numel(mode)<1  % nothing to insert
    return
end;

modes={'AND' 'OR' 'OVER' 'OFF' ''};
if ~any(strcmp(mode,modes))
    error(['Invalid mode string: ' mode]);
end;
mask.merge=mode;
mask.encoding='RLE';  % RLE
mask.data=RunLengthEncode(msk);

if ~isfield(mi,'mask') || numel(mi.mask)<1
    mi.mask=struct('merge',[],'encoding',[],'data',[]);
end;
mi.mask(index)=mask;
