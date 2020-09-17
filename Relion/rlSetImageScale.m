function [mi,m,origImageSize]=rlSetImageScale(mi,mode,nFrames);
% Set the size and scale parameters in the mi struct, possibly loading and
% scaling the image m as well.
% mi.imageSize is always set to the "nice" padded image size. If you are
% going to use the original micrograph as the working image, you should 
% change this to mi.imageSize=origImageSize.
% There are three modes.
% 1. Read the micrograph, pad the image to a nice size, convert to fractional contrast
%   (this is our traditional method.) Set mi.imageSize and mi.imageMedian.
%   Scale up by sqrt(nFrames), which should be >1 to fix MotionCor2's wrong
%   dc scaling.
% 2. Read the micrograph, insert its actual size as mi.imageSize, and
%   compute mi.imageNormScale and mi.imageMedian.
% 3. Estimate mi.imageNormScale from cpe,pixA, and mi.doses(1). 
%     In modes 1 and 2, the returned image m is padded to the nice size and
%     already normalized.
% For mode 3:
%     Return m=0; Assume that we'll use the original micrograph as the
%     "processed" or "merged" image. Compute the normalization we'll have to
%     apply to make the AC image component approximately equal to fractional
%     contrast. We'll estimate the variance by averaging the 1D power spectrum
%     from 0.3 to 0.7 * Nyquist. Assuming an arbitrary image scaling by a, the
%     est variance of a counting image should be a^2 * pixelDose, which will be
%     (a * size of original pixel)^2 * doses(1), with a being the unknown
%     scaling factor. (By the scaling used by MotionCor2, the "size of original
%     pixel" should be a superres pixel. The final fractional-contrast image
%     will have the variance 1/pixelDose. We'll have to get it by scaling the
%     raw image m0 by 1/sqrt(est variance*pixelDose). In the end we'll get the
%     scaled micrograph, after computing the median, by
%     scaledImg = ( m0-median(m0(:)) )*mi.imageNormScale;

    m0=0; % default returned value.
    m=0;
micName=[mi.imagePath mi.imageFilenames{1}];
micFound=exist(micName,'file');
if mode<3
    if micFound
        m0=ReadMRC(micName);
        med=median(m0(:));
    else
        disp(['Micrograph file not found: ' micName]);
        mode=3;
        med=0;
    end;
end;
origImageSize=size(m0);
if numel(origImageSize)<2
    origImageSize(2)=origImageSize(1);
end;
niceImageSize=NextNiceNumber(origImageSize,5,8);
mi.imageSize=niceImageSize;
mc=Crop(m0-med,niceImageSize);

switch mode
    case 1 % We assume we know the DC component correctly, or else use the
        % correction for MotionCor2's error to scale up [should it be down?]
        % the DC value.
        mi.imageNormScale=1/med*sqrt(nFrames);
        mi.imageMedian=med;
        m1=RemoveOutliers(mc);
        m=m1*mi.imageNormScale;
    case 2 % don't know scaling or the DC value, estimate from the image.
%       Estimate shot noise from the mean power spectrum 0.3 ... 0.7 x Nyquist
        sds=floor(min(mi.imageSize)/256); % Downsampling factor for spectrum
        mCtr=single(Crop(m0,256*sds)); % Grab a square region of the image.
        sp1=RadialPowerSpectrum(mCtr,0,sds);
        spN=numel(sp1);
        spLims=round(spN*[.3 .7]);
        estVar=sum(sp1(spLims(1):spLims(2)))/diff(spLims);
        mi.imageNormScale=1/(mi.pixA*sqrt(mi.doses(1)*estVar));
%         Multiply the image by this to get a normalized image.
        mi.imageMedian=med; % best estimate we have of the main image mean.
        m=mc*mi.imageNormScale;

    case 3 % Don't have an image available
    %      Image statistics: in the absence of image data,
    %         we assume a true counting camera, binned from superres, and
    %         just wing it for scale. Don't return m.
        mi.imageNormScale=1/(mi.cpe*mi.pixA^2*mi.doses(1));
        % mi.med=0; % don't assign image median.
    otherwise
        error('Mode not recognized');
end;