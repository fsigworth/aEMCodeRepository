function mc=Cropz(m,n,isstack,fillval)
% function mc=Cropz(m,n,isstack,fillvalue)
% Reduce the size of the 1d array, square or cube m by cropping
% (or increase the size by padding with fillval, by default zero)
% to a final size of n x n or n x n x n by removing only the outermost points.
% This is similar to Crop() but does not affect points near (1,1).
% Thus m1=Cropz(m,[nx ny]) give the same result as
% m1=m(1:nx,1:ny), with any empty space filled with fillval.
% If m is 2-dimensional and n is a vector, m is cropped to n=[nx ny].
% If the flag isstack = 1 then a 3D array m is treated as a stack of 2D
% images, and each image is cropped to n x n.
% For 2D images, the input image doesn't have to be square.
% The result is double if fillval is double; by default the result is
% single.
% Now handles simultaneous cropping/padding for 2D input. fs Apr-20.
% 3D-4D inputs are still
% assumed to be stacks of square, or cubic in dimension.

if nargin<3
    isstack=0;
end;
if nargin<4
    fillval=single(0);  % Force a single output when padding.
else
    fillval=single(fillval);
end;

sz=size(m);
ndi=ndims(m);
if ndi==2 && any(sz==1)
    ndi=1;
end;

switch ndi
    case 1
        n1=numel(m);
        n=max(n); % force a column vector of output.
        m=reshape(m,n1,1); % force a column vector of input
        ns=n1-n;  % positive for cropping, neg for padding.
        if ns<0
            mc=fillval*ones(n(1),1); % padding
            mc(1:n1(1),1)=m;
        else
            mc=m(1:n1(1),1);
        end;
         
    case 2
        if numel(n)<2
            n(2)=n(1);
        end;
        no=min(sz,n);
        if any(n>sz)
            mc=fillval*ones(n);
        end;
        mc(1:no(1),1:no(2))=m(1:no(1),1:no(2));
        
    case 3 % m is 3D
        if isstack % a stack of 2D images
            if numel(n)<2
                n=n(1)*[1 1];
            end;
            ns=sz(1:2)-n(1:2);  % Shift term for scaling down.
            n(3)=sz(3);
            if all(ns>=0) % cropping down
                mc=m(1:ns(1)+n(1),1:ns(2)+n(2),:);
            elseif all(ns<=0)  % padding
                mc=fillval*ones(n,'single');
                mc(1:sz(1),1:sz(2),:)=m;
            else
                error('Can''t crop and pad the same image');
            end;
        else  % not a stack
            if numel(n)<3
                n=n(1)*[1 1 1];
            end;
            ns=floor(sz/2)-floor(n/2);  % Shift term for scaling down.
            if all(ns>=0) % cropping down
                mc=m(1:n(1),1:n(2),1:n(3));
            elseif all(ns<=0)
                mc=fillval*ones(n);
                mc(1:sz(1),1:sz(2),1:sz(3))=m;
            else
                error('Can''t crop and pad dimensions simultaneously');
            end;
        end;
    case 4
        % this must be a stack of 3D
            if numel(n)<3
                n=n(1)*[1 1 1];
            end;
            n(4)=sz(4);
            if all(sz>=n) % cropping down
                mc=m(1:n(1),1:n(2),1:n(3),:);
            elseif all(sz<=n)
                mc=fillval*ones(n);
                mc(1:sz(1),1:sz(2),1:sz(3),:)=m;
            else
                error('Can''t crop and pad dimensions simultaneously');
            end;
    otherwise
        error(['Cropz: dimension too large: ' num2str(ndi)]);
end;
