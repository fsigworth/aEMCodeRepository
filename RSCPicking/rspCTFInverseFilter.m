function m1=rspCTFInverseFilter(m,mi,compFraction,pars)
% function m=rspCTFInverseFilter(m,mi,compFraction,pars)
% Given the downsampled image m, boost the low frequencies below the 1st
% peak of the first exposure's CTF.  If compFraction=1, the boost is the inverse of
% the ctf, up to a value of 10.
% If we're padding the image m, it would be good to have it correspond to
% mi.padImageSize.
% Includes code to make a fake double-exposure image if weights are not all 1.
%
% mi fields that are used
%   weights (if doDoubleExposure only)
%   pixA
%   padImageSize (or imageSize if padImageSize not given)
%   ctf

if nargin<4
    pars=struct;
end;
% Here are default parameters
dPars.phaseFlip=1;
dPars.firstPeakAmp=.5;  % Constrained amplitude during first peak
dPars.k=.05;
dPars.fHP=0; % Nonzero to give a highpass filter too.
dPars.ds=0; % If nonzero, forces the downsampling factor
dPars.doDoubleExposure=0; % Emulate double exposure using only the first CTFf

if nargin<4
    pars=struct;
end;

pars=SetDefaultValues(dPars,pars,1);

if compFraction<.01 && ~pars.doDoubleExposure % nothing to do
    m1=m;
    return
end;

n=size(m);

% Figure out the downsampling factor
if pars.ds==0
    if ~isfield(mi,'padImageSize');
        origSize=mi.padImageSize;
    else
        origSize=mi.imageSize;
    end;
    ds=origSize(1)/n(1);  % actual downsampling of m
end;

pixA=mi.pixA*ds;

freqs=RadiusNorm(n)/pixA;  % frequencies of padded image
mi.ctf(1).B=0; % ignore envelope function
[ctf,chi]=ContrastTransfer(freqs,mi.ctf(1));  % Get Chi too.
maxChi=1-asin(pars.firstPeakAmp)/pi;
broadChi=1-asin(0.01)/pi; % all the way up to the first zero
peakMask=chi>-maxChi;  % mask up to just beyond the first peak
broadMask=chi>-broadChi;
% mxAmp=max2d(peakMask.*abs(ctf));
H=ones(n,'single');

%     Make H such that c.*H is ~flat up to maxChi, no boost beyond.
    H(peakMask)=(pars.firstPeakAmp+pars.k)./(abs(ctf(peakMask))+pars.k);
    H(broadMask)=(1-compFraction)+compFraction*H(broadMask);

if pars.phaseFlip
    H=H.*sign(-ctf);
end;
H=H.*GaussHPKernel(n,pars.fHP); % does nothing if fHP==0

if pars.doDoubleExposure
    c0=meGetEffectiveCTF(mi,n,ds);
    mi.weights=effWeights;
    c1=meGetEffectiveCTF(mi,n,ds);
    H=H.*c1./c0;
%    H=H.*(1+weightComp*(c1./c0-1));
% figure(2); clf; plot([sectr(c1) sectr(c0) sectr(c1./c0)]); figure(1);
end;
%%%%%%%%%

m1=real(ifftn(fftn(m).*ifftshift(H)));
