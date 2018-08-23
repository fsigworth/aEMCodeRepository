function modelParticle=rspMakeModelParticle(dis,mi,coords,rscc)
% Based on cross-correlation values at the coords location in rscc, create
% a model particle from eigenimages.  The returned modelParticle image is
% the same size as an eigenimage.

[nd,~,nt]=size(rscc.eigenImgs);
ix=round(coords(1)/dis.ds+1); % Coordinates in our images
iy=round(coords(2)/dis.ds+1);
templ=single(rscc.mxTemplInds(ix,iy)); % template number
if templ<1
    modelParticle=zeros(nd,nd,'single');
    return
end;

ampU=rscc.mxCCU(ix,iy);  % get the absolute amplitude
if ampU==0 % we use the relative amplitude instead
%     disp('Using scaled amplitude');
    iv=max(1,rscc.mxVesInds(ix,iy));
    ampU=mi.vesicle.s(iv)*rscc.mxCC(ix,iy);
end;
% Perform the matrix multiplication eigenImages*vList
vsz=numel(rscc.vList)/nt;
vlst=reshape(rscc.vList,nt,vsz);
eigs=reshape(rscc.eigenImgs,nd^2,nt);
modelParticle=-ampU*reshape(eigs*vlst(:,templ),nd,nd);
