function moix=reRescaleMoi(moi,ri)
% function moix=reRescaleMoi(moi,ri,gFlags)
% Down- or up-sample the following fields of the moi or gMoi structure:
%   refVols
%   ringFits (only those flagged by gFlags if given)
%   activeTrans (nearest-neighbor, sized according to maxShift)
% This is used in changing resolution of reconstructions.
% We use ri.nCurrent and ri.nTrans fields from ri, and assume that ri
% reflects the new size.

% Rescaling of sigmaN.  Theoretically should be 1, but this is conservative
noiseExp=1.5;

n=ri.nCurrent;
nt=ri.nTrans;

n0=size(moi.refVols,1); % original image size
us=n/n0;  % upsampling ratio

totalImgs=size(moi.imgAmps,1);

nv=size(moi.refVols,4);
moix=moi;
moix.refVols=zeros(n,n,n,nv,'single');
for i=1:nv
    moix.refVols(:,:,:,i)=DownsampleGeneral(moi.refVols(:,:,:,i),n);
end;

if ri.flags.mode2D
    moix.refs=DownsampleGeneral(moi.refs,n,[],1);
end;

% Resample the active translations
nt0=sqrt(size(moi.activeTrans,1));
moi.activeTrans=reshape(moi.activeTrans,nt0,nt0,totalImgs);
moix.activeTrans=DownsampleNearest(moi.activeTrans,nt,0,1);
moix.activeTrans=reshape(moix.activeTrans,nt^2,totalImgs);
disp(['Upsampling moi.activeTrans from ' num2str(nt0) ' to ' num2str(nt)]);


% Resample the active alphas
ni=size(moi.activeAlphas,2);
na0=size(moi.activeAlphas,1);  % original nAlphasI
na=size(ri.alphas,1);  % new number of alphas
moix.activeAlphas=false(na,ni,totalImgs);

actA=moi.activeAlphas;

disp(['Upsampling moi.activeAlphas from ' num2str(na0) ' to ' num2str(na)]);
for i=1:totalImgs
    for j=1:ni
        a=DownsampleNearest(actA(:,j,i),na);
        a(1:end-1)=a(1:end-1) | a(2:end);
        a(2:end)=a(2:end) | a(1:end-1);
        moix.activeAlphas(:,j,i)=a;
    end;
end;

if ~ri.flags.mode2D % for 3D reconstruction, resample gamma and beta
    sz=size(moi.activeRVs);
    ngb=ri.angleN(3:-1:2);  % gamma,beta
    moi.activeRVs=reshape(moi.activeRVs,[sz(1:2) prod(sz(3:end))]);
    moix.activeRVs=DownsampleNearest(moi.activeRVs,ngb,0,1);
    moix.activeRVs=BinaryConvolve(moix.activeRVs,true(3,3));
    moix.activeRVs=reshape(moix.activeRVs,[ngb sz(3:end)]);
    disp(['Upsampling moi.activeRVs from ' num2str(round(sum(moi.activeRVs(:))/totalImgs))...
        ' to ' num2str(round(sum(moix.activeRVs(:))/totalImgs))]);
end;

% Change the intensive fields
moix.sigmaN=moi.sigmaN*(us.^noiseExp);
moix.sigmaC=moi.sigmaC*us;
moix.sigmaG=moi.sigmaG*us;
% moix.b0=moi.b0*us;

