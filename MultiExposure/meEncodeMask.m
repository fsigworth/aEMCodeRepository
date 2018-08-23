function mask=meEncodeMask(img,mode)
% Given a binary image msk (masked points are zeros), encode it into a
% structure that can be copied e.g.
% mi.mask(i)=mstruct;
% This information is decoded by meGetMask.

if nargin<2
    mode='AND';
end;

modes={'AND' 'OR' 'OVER' 'OFF' ''};
if ~any(strcmp(mode,modes))
    error(['Invalid mode string: ' mode]);
end;

mask.merge=mode;
mask.encoding='RLE';  % RLE
mask.data=RunLengthEncode(img);

% if ~isfield(mi,'mask')
%     mi.mask=struct('merge',[],'encoding',[],'data',[]);
% end;
% mi.mask(index)=mask;
