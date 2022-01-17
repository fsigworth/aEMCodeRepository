function [H,filterOk]=meGetNoiseWhiteningFilter(mi,n,ds,nzeros,fHP,effCTF)
% function [H,filterOk]=meGetNoiseWhiteningFilter(mi,n,ds,nzeros,fHP,effCTF)
% Given the micrograph info structure mi, returns a filter frequency
% response H of size n, f=0 in the center, that is the inverse of the
% mi.noiseModel as filtered by the effective CTF. 
% The filter is scaled such that the HF asymptote is unity:
%     H=sqrt(mean(shot(:))./(spec.*effCTF.^2+shot)); % =1 at high frequencies
%
% If fHP is given and >0,
% an additional Gaussian HP with half-power frequency fHP (in A^-1) is applied.
% The result is for an image downsampled by ds (default ds=mi.imageSize/n).
% No need to give ds, or ds=0 is ok, if the image wasn't cropped.
% If effCTF is not given, it is computed from meComputeMergeCoeffs2.
%    H=sqrt(shot/(excessSpectrum * ctf^2 + shotSpectrum)).
% filterOk is true if the mi structure contained valid filter parameters.
if nargin<3
    ds=0;
end;
if nargin<4
    nzeros=1;
end;
if nargin<5
    fHP=0;
end;
if ~isfield(mi,'weights')
    mi.weights=mi.doses>0;
end;

if ds==0
    ds=mi.imageSize(1)/n(1);
end;

pixA=mi.pixA*ds;  % Downsampled pixel size
f2d=RadiusNorm(n)/pixA;  % frequency in A^-1, zero center

filterOk=isfield(mi,'noiseModelPars') && numel(mi.noiseModelPars)>0 && ~any(isnan(mi.noiseModelPars));
if filterOk
    % Get the effective CTF and noise spectrum
    if nargin<6 % we need to get the effCTF
        [~, effCTF]=meComputeMergeCoeffs2(f2d,mi,nzeros);
    end;
    [spec, shot]=meEvalNoiseModel(f2d,mi);
    spec=max(spec,0);  % This excess noise should be positive
    H=sqrt(mean(shot(:))./(spec.*effCTF.^2+shot)); % =1 at high frequencies
else
    H=1;
end;
if fHP>0
    k=(log(2)/2)*fHP.^2;
    H=H.*exp(-k./(f2d.^2+1e-9));	% inverse Gaussian kernel
end;
H(isnan(H))=0;  % zero out any NaN values.