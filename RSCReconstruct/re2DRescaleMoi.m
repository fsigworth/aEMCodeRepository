function moix=re2DRescaleMoi(moi,ri)
% function moix=reRescaleMoi2D(moi,ri)
% Version for 2D classification
% Down- or up-sample the following fields of the moi or gMoi structure:
%   refCls
%   ringFits (only those flagged by gFlags if given)
%   activeTrans (nearest-neighbor, sized according to maxShift)
% This is used in changing resolution of reconstructions.
% We use ri.nCurrent and ri.nTrans fields from ri.
n=ri.nCurrent;
nt=ri.nTrans;

n0=size(moi.refs,1); % original image size
us=n/n0;  % upsampling ratio

totalImgs=numel(moi.imgAmps);

nRefs=size(moi.refs,3);
moix=moi;
moix.refs=zeros(n,n,nRefs,'single');
for i=1:nRefs
    moix.refs(:,:,i)=DownsampleGeneral(moi.refs(:,:,i),n);
end;

% Resample the active translations
nt0=sqrt(size(moi.activeTrans,1));
moi.activeTrans=reshape(moi.activeTrans,nt0,nt0,totalImgs);
moix.activeTrans=DownsampleNearest(moi.activeTrans,nt,0,1);
moix.activeTrans=reshape(moix.activeTrans,nt^2,totalImgs);

% Change the intensive fields
moix.sigmaN=moi.sigmaN*us;
moix.sigmaC=moi.sigmaC*us;
moix.sigmaG=moi.sigmaG*us;
% moix.b0=moi.b0*us;

