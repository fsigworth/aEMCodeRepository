function [mxVals,mxInds,eigCCs]=hrCCSearch(mc,ei);

fM=fftn(ifftshift(mc));
n1=size(mc);
n=size(ei.eigImgs,1);
nv=size(ei.eigImgs,3);
nRefs=size(ei.normCoeffs,1);
eigCCs=zeros([prod(n1) nv],'single');

disp([num2str(nv) ' cross correlations']);
for i=1:nv
    fRef=conj(fftn(Crop(ei.eigImgs(:,:,i),n1)));
    eigCCs(:,i)=reshape(real(ifftn(fRef.*fM)),prod(n1),1);
    if mod(i,100)==0
        fprintf('.'); % a dot every 100 eigenimages
    end;
end;
fprintf('\n');

%%
disp(['Expansion to ' num2str(nRefs)]);
% at nv=503, this runs about 4x faster than brute-force fft calc.
mxV=-ones(prod(n1),1,'single')*inf;
mxI=zeros(prod(n1),1,'int32');
nBlock=100;
ind=0;
while ind<nRefs
    block=min(nBlock,nRefs-ind);
    blockCCs=eigCCs*ei.normCoeffs(ind+1:ind+block,:)'; % n1^2 x nRefs
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

mxVals=reshape(mxV,n1);
mxInds=reshape(mxI,n1);
eigCCs=reshape(eigCCs,[n1 nv]);
