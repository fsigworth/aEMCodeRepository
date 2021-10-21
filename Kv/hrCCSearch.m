function [mxVals,mxInds,eigCCs]=hrCCSearch(mc,ei,upFactor)
if nargin < 3
    upFactor=1; % upsampling factor
end;
fM=fftn(ifftshift(mc));
n1=size(mc);
nw=n1*upFactor;
n=size(ei.eigImgs,1);
nv=size(ei.eigImgs,3);
nRefs=size(ei.normCoeffs,1);
eigCCs=zeros([prod(n1) nv],'single');

disp([num2str(nv) ' cross correlations']);
for i=1:nv
    fRef=conj(fftn(Crop(ei.eigImgs(:,:,i),n1)));
%     if upFactor>1
%         cc=real(ifftn(Crop(fRef.*fM,nw))); % oversampling
%     else
        cc=real(ifftn(fRef.*fM));
%     end;
    eigCCs(:,i)=reshape(cc,prod(n1),1);
    if mod(i,100)==0
        fprintf('.'); % a dot every 100 eigenimages
    end;
end;
fprintf('\n');

%%
disp(['Expansion to ' num2str(nRefs)]);
% at nv=503, this runs about 4x faster than brute-force fft calc.
mxV=-ones(prod(nw),1,'single')*inf;
mxI=zeros(prod(nw),1,'int32');
nBlock=NextNiceNumber(100/upFactor^2,7,-1);
ind=0;
while ind<nRefs
    block=min(nBlock,nRefs-ind);
    blockCCs=eigCCs*ei.normCoeffs(ind+1:ind+block,:)'; % nw^2 x nRef
    if upFactor>1  % upsample the resulting block CCs
        blockCCs=reshape(Downsample(reshape(blockCCs,[n1 block]),nw,1),prod(nw),block);
    end;
    [bMxV,bIndV]=max(blockCCs,[],2);
    blockBigger=bMxV>mxV;
    mxV(blockBigger)=bMxV(blockBigger);
    mxI(blockBigger)=bIndV(blockBigger);
    ind=ind+block;
    if mod(ind,1000)==0
        fprintf('.'); % a dot every 1000 eigenimages.
    end;
end;
fprintf('\n');

mxVals=reshape(mxV,nw);
mxInds=reshape(mxI,nw);
eigCCs=reshape(eigCCs,[nw nv]);
