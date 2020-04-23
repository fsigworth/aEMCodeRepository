function msk=meGetMask(mi,n,indices)
% function msk=meGetMask(mi,n,indices)
% From the mi.mask array of structures, generate a boolean image of size n.
% If the optional argument indices is given, use only those elements of the
% mask stack.
% The masked points have the value logical false.
% If no masks are present, return the default msk = true(n).
% If n=0 then we return the mask of index indices(1) at its native size.
% (we can return only one cuz we can't guarantee the various masks are
% defined at the same sizes).
% 
% mi.mask has the fields
%   merge  - a text string AND, OR, OVER, OFF
%   encoding - a text string RLE, DISC
%       where RLE is run-length encoding of a binary image, and RIM is
%       the exterior of a binary disc.
%   data - data encoding the mask.
%   for 'RIM', data = [X Y R], all in normalized units; X,Y=0 is in the
%   center of the image; R=0.5 just touches the edges.
%   for 'RLE' data is a structure that is read by the RLE decoder.

maxExpandedImage=3*1024;  % maximum temporary array dimension;
%                       if larger than this, we use Fourier interpolation
if n==0
    if nargin<3 || numel(indices)>0
        indices=indices(1);
    else
        error('No index value given for single mask');
    end;
    msk=true;
else
    msk=true(n);  % default
end;

if isfield(mi,'mask')
    nim=numel(mi.mask);
    if nargin<3 || numel(indices)<1
        indices=1:nim;
    else
        if ~any(indices<=nim)
            return
        end;
        indices=indices(indices<=nim);
    end;
    for i=indices
        ms=mi.mask(i);
        if numel(ms.merge)>0  && numel(ms.encoding)>0 ...
                        && ~strcmpi(ms.encoding,'OFF')  % something there
            switch mi.mask(i).encoding
                case 'RLE'
                    m1=RunLengthDecode(mi.mask(i).data);
                case {'RIM' 'beam'}
                    ctr=ceil((n+1)/2);
                    d=mi.mask(i).data;
                    m1=fuzzymask(n,2,d(3)*n(1),0,d(1:2).*n+ctr);
                otherwise
                    error(['Unexpected mi.mask.encoding: ' mi.mask(i).encoding]);
            end;
            if n==0
                msk=m1;
                return
            end;

            n1=size(m1);
            if n1(1) ~= n(1)  % need to change size
                lcf=LeastCommonFactor(n(1),n1(1));
                if n1(1)*n(1)/lcf <= maxExpandedImage  % don't allow huge expansions
                    m1=BinImage(ExpandImage(m1,n(1)/lcf),n1(1)/lcf);
                else
                    m1=DownsampleGeneral(m1,n,n(1)/n1(1))>0.5;
                end;
            end;
            if any(size(m1)~=n)
                disp('Mismatched mask dimensions; cropping.');
                m1=Crop(m1,n);
            end;
            switch mi.mask(i).merge
                case 'AND'
                    msk=msk & m1;
                case 'OR'
                    msk=msk | m1;
                case 'OVER' % Overwrite the earlier masks.
                    msk=m1;
                case 'OFF'
                    % do nothing
                otherwise
                    error(['Unexpected mi.mask.merge value: ' mi.mask(i).merge]);
            end;
        end;
    end;
else
    msk=true(n);  % default is no mask.
end;