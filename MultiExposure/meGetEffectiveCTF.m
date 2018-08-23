function c=meGetEffectiveCTF(mi,n,ds,mergeMode,nzeros)
% function c=meGetEffectiveCTF(mi,n,ds,mergeMode,nzeros)
% Given the micrograph info structure mi, the size n of the desired 2D ctf
% and the downsampling ratio ds, return the effective ctf including the
% effects of merging and CCD (after doing CCD pre-whitening).
% By default, ds=mi.imageSize(1)/n.  You shouldn't use the default if n
% represents a cropped image, e.g. a particle image. The default is used if
% ds is given as 0.  The returned c has zero frequency in the center.  By
% default c is non-negative.  mergeMode selects the various merging modes in meComputeMergeCoeffs2().
% if mergeMode=3 there is no phase-flipping, and c is based on the first image's raw ctf.
% fs 3 Sep 11, 15 Aug 16

if numel(n)<2
    n=[1 1]*n;
end;
if nargin<3 || ds==0
    ds=mi.imageSize(1)/n(1);
end;
if nargin<4
    mergeMode=1;
end;
mergeMode=max(1,mergeMode);

if nargin<5
    nzeros=1;  % default is to include only the 1st zero in later images
end;
if ~isfield(mi,'weights') || numel(mi.weights)<numel(mi.doses)
    mi.weights=single(mi.doses>0);
end;
effPixA=mi.pixA*ds;
freqs=RadiusNorm(n)/effPixA;  % frequencies for evaluating the CTF
% Change the weights to include the ampFactors
weights=mi.weights;
if isfield(mi.ctf,'ampFactor') && numel(mi.ctf)>=numel(weights)
    for i=1:numel(weights)
        weights(i)=weights(i)*mi.ctf(i).ampFactor;
    end;
end;
[~, effctf]=meComputeMergeCoeffs3(freqs, mi,nzeros,mergeMode);
% if noflip
%     effctf=effctf.*sign(coeffs(:,:,1));
% end;
c=effctf.*CCDEffCTF(mi.camera,n,ds);
