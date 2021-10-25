function [mxVals,mxInds]=hrProjSearch(mc,projs,upFactor)
if nargin < 3
    upFactor=1; % upsampling factor
end;
n1=size(mc);
nw=n1*upFactor;
nRefs=size(projs,3);
fM=fftn(ifftshift(mc));

mxV=-ones(prod(nw),1,'single')*inf;
mxI=zeros(prod(nw),1,'int32');
nBlock=NextNiceNumber(100/upFactor^2,7,-1);
blockCCs=zeros([nw nBlock],'single');
ind=0;
first=1;
tic;
while ind<nRefs
    if nRefs-ind<nBlock
        nBlock=nRefs-ind;
        blockCCs=zeros([nw nBlock],'single');
    end;
%     block=min(nBlock,nRefs-ind);
    for j=1:nBlock
        fRef=conj(fftn(Crop(projs(:,:,ind+j),n1)));
        blockCCs(:,:,j)=real(ifftn(Cropo(fRef.*fM,nw,1)));
    end;
    [bMxV,bIndV]=max(reshape(blockCCs,[prod(nw) nBlock]),[],2);
    blockBigger=bMxV>mxV;
    mxV(blockBigger)=bMxV(blockBigger);
    mxI(blockBigger)=bIndV(blockBigger);
    if mod(ind,1000)<nBlock
        if first && ind>=1000  % print out time for 1000 ccs
            fprintf(['For ' num2str(ind) ' ccs ']);
            toc;
            first=0;
        else
        fprintf('.'); % a dot every 1000 eigenimages.
        end;
    end;
    ind=ind+nBlock;
end;
fprintf('\n');

mxVals=reshape(mxV,nw);
mxInds=reshape(mxI,nw);
