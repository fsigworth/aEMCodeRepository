function [mi,m]=rlSetImageScale(mi,mode,nFrames);
% Set the size and scale parameters in the mi struct, possibly loading and
% scaling the image m as well.
% There are three modes.
% 1. Read the micrograph, pad the image to a nice size, convert to fractional contrast
%   (this is our traditional method.) Set mi.imageSize and mi.imageMedian
% 2. Read the micrograph, insert its actual size as mi.imageSize, and
%   compute mi.imageNormScale and mi.imageMedian.
% 3. Estimate mi.imageNormScale from cpe,pixA, and mi.doses(1).
% For modes 2 and 3:
%     Set up to use the original micrograph as the "processed" or "merged"
%     image.
%     Compute the normalization we'll have to apply to make the AC image component
%     approximately equal to fractional contrast. We'll estimate the variance
%     by averaging the 1D power spectrum from 0.3 to 0.7 * Nyquist. Assuming
%     an arbitrary image scaling by a, the est variance of a counting image should be
%     a^2 * pixelDose, which will be (a * size of original pixel)^2 * doses(1),
%     with a being the unknown scaling factor. (By
%     the scaling used by MotionCor2, the "size of original pixel" should be a
%     superres pixel. The final fractional-contrast image will
%     have the variance 1/pixelDose. We get it by scaling the raw image by
%     1/sqrt(est variance*pixelDose).
%     In the end we'll get the scaled micrograph as
%     scaledImg = (m-mi.imageMedian)*mi.imageNormScale;

    m=0; % default returned value.

micName=[mi.imagePath mi.imageFilenames{1}];
micFound=exist(micName,'file');
if mode<3
    if micFound
        m0=ReadMRC(micName);
        med=median(m0(:));
    else
        disp(['Micrograph file not found: ' micName]);
        mode=3;
    end;
end;

switch mode
    case 1
        nsz=NextNiceNumber(size(m0),5,8);
        m=Crop((m0-med)/med*sqrt(nFrames),nsz);
        mi.imageSize=nsz;
    case 2
        mi.imageSize=size(m0);
        sds=floor(min(mi.imageSize)/256); % Downsampling factor for spectrum
        mc=single(Crop(m0,256*sds)); % Grab a square region of the image.
        sp1=RadialPowerSpectrum(mc,0,sds);
        spN=numel(sp1);
        spLims=round(spN*[.3 .7]);
        estVar=sum(sp1(spLims(1):spLims(2)))/diff(spLims);
        mi.imageNormScale=1/(mi.pixA*sqrt(mi.doses(1)*estVar));
%         Multiply the image by this to get a normalized image.
        mi.imageMedian=med; % best estimate we have of the main image mean.
        m=(single(m0)-med)*mi.imageNormScale;
    case 3
    %      Image statistics: in the absence of image data,
    %         we assume a true counting camera, binned from superres, and
    %         just wing it for scale. Don't return m.
        mi.imageNormScale=1/(mi.cpe*mi.pixA^2*mi.doses(1));
    otherwise
        error('Mode not recognized');
end;