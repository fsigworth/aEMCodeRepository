function [out,vars]=NormalizeImages(m,varNormalize,applyMask,globalNorm)
% function out=NormalizeImages(m,varNormalize,applyMask)
% Normalize each image in an image stack for classification and alignment.
% Images are assumed to be square.
% Defaults are varNormalize=1 and applyMask=1.
% We use an annulus for estimating variance and mean value.  We normalize
% it so that the annulus mean is zero and variance is unity.  If mask=1 a
% circular mask is also applied, 95% of max radius.
if nargin<4
    globalNorm=0;
end;
if nargin<3
    applyMask=1;
end;
if nargin<2
    varNormalize=1;
end;
OutR=0.5;  % Outer radius
InR=0.4;   % Inner radius
Rise=0.05;   % Feathering of masks

[n ny nim]=size(m);

% Make the masks for normalizing images
Masko=fuzzymask(n,2,OutR*n,Rise*n); % Outer mask at 90% of radius
Maski=fuzzymask(n,2,InR*n,Rise*n); % Inner mask at 75% of radius
AnnMask=Masko-Maski;  % Annular mask for measuring mean and variance
AnnMaskSum=sum(AnnMask(:));

out=single(zeros(n,n,nim));
vars=zeros(nim,1);
am2=AnnMask(:)'*AnnMask(:);
for i=1:nim
    m0=double(m(:,:,i));
    m1=AnnMask.*m0;  % Pick out the annulus to estimate mean and variance
    mbkg=sum(m1(:))/AnnMaskSum;
    msub=m0-mbkg;
    m2=AnnMask.*msub;
    var=m2(:)'*m2(:)/am2;
    
    % Subtract the constant background, apply outer mask
    if applyMask
        m1=Masko.*msub;
    else
        m1=msub;
    end;
    if ~globalNorm && varNormalize && var>1e-12 % avoid division by zero
        m1=m1/sqrt(var);
    end;
    vars(i)=var;
    out(:,:,i)=m1;
end;
if globalNorm
    out=out/sqrt(median(vars));
end;
