function ccs=spRotCorr(pims,pref,doFlip)
% Given polar-transformed reference pref and images pims, obtain the
% rotational cross-correlations (image relative to ref)
% The cross correlations are in the first dimension of the 
% ntheta x nims x nrefs x nflips returned variable ccs.
%  Example of quick rotational alignment
    % img=GaussFilt(randn(64,64),.2);
    % imgs=single(zeros(64,64,nim));
    % refs=single(zeros(64,64,nim));
    % for i=1:nims
    %     imgs(:,:,i)=grotate(img,(i-1)*0.1)+5*randn(64,64);
    %     refs(:,:,i)=grotate(img,(i-1)*.01)+randn(64,64);
    % end;
    % pimgs=gridToPolar(imgs);
    % prefs=gridToPolar(refs);
    % ccs=spRotCorr(pimgs,prefs);



if nargin<3
    doFlip=1;
end;

[nr nt nims]=size(pims);
[nrr ntr nrefs]=size(pref);
weights=(0.5:1:nr)'/nr;  % weight of each pixel
nflip=1+(doFlip>0);
ccs=zeros(nt,nims,nrefs);

fims=fft(pims,[],2);  % take the ffts along theta
fref=conj(fft(pref,[],2));
if nflip>1
    fref(:,:,:,2)=conj(fref);
end;
for i=1:nims
    for j=1:nrefs
        for k=1:nflip
%             ccs(:,i,j,k)=weights*real(ifft(fims(:,:,i).*fref(:,:,j,k),[],2));
   ccs(:,i,j,k)=real(ifft((fims(:,:,i).*fref(:,:,j,k))'))*weights;
        end;
    end;

end;


