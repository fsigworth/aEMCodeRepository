function out=ExpandImage(in,nb)
% function out=ExpandImage(in,nb)
% Reverse of BinImage.  Presently works with 1D vectors and with 2D images or stacks.
% The class of out is the same as in.

if nb<=1
    out=in;
    return;
end;
nb=round(nb);

dims=ndims(in);
[nx, ny, nim]=size(in);
n=[nx ny];
dims=2;
if any(n==1)  % a 1d function
    dims=1;
    n=prod(n);
    %     in=reshape(in,n,nim);
end;
switch dims
    case 1
        in1=repmat(in(:)',nb,1);
        out=in1(:);
    case 2
        if isnumeric(in)
            out=zeros([n.*nb nim],'like',in);
        elseif islogical(in)
            out=false([n.*nb nim]);
        else
            error('Input is not numeric or logical');
        end;
            
        for i=1:nim
            if numel(nb)<2  % if nb is a scalar, make it a vector
                nb=nb(1)*ones(1,2);
            end;
            
            in1=in(:,:,i)'; % First, take the transpose
            q=repmat(in1(:)',nb(2),1);  % expand in y direction
            qr=reshape(q,ny*nb(2),nx)'; % transpose back
            out1=repmat(qr(:)',nb(1),1); % expand in x direction
            out(:,:,i)=reshape(out1,nx*nb(1),ny*nb(2));
        end;
    otherwise
        error('Only 2d images handled.');
end;
