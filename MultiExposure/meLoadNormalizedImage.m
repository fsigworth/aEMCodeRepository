function [mOut,M,ok,mIn]=meLoadNormalizedImage(mi,target,sufExts)
% function [mOut,M,ok,mIn]=meLoadNormalizedImage(mi,targetSize,sufExts)
% function [mOut,M,ok,mIn]=meLoadNormalizedImage(mi,maxPixA,sufExts)
%  -targetSize is a two-element vector.
%  -maxPixA is a scalar.
%  -M is the affine transform of mOut's possible downsampling and shift,
% Loads one of the variant merged image files, or the original micrograph
% in the case that we have normalization information in the mi structure.
% Return an image mOut that is cropped and downsampled to fit within
% targetSize. Alternatively, a nice number output size is chosen 
% to have a pixel size <= maxPixA and an integer downsampling factor.
% The returned matrix M has the elements
%   M(1,1): downsampling factor from original micrograph to mOut
%   -M(1:2,3): shift added to original micrograph before downsampling. Thus to
%   determine a pixel coordinate in the the original micrograph, it is
%   computed so:
% origMicrographXY=M*[xOut;yOut;1];
%  If desired, the full-size imported image is available as mIn.

if nargin<3
    sufExts={'s.mrc' 'z.tif' '.mrc'}; % default: prefer small or compressed files.
end;
ok = false;
M=zeros(3,3);

% Load the merged image. We assume it's already scaled to unit intensity but
% with the mean (background) value subtracted. This means for image padding
% we can use zeros.
imageBasename=[mi.procPath mi.baseFilename 'm.mrc'];
[fullMergedImageName,nameOk]=CheckForAltImage(imageBasename,sufExts);  % valid filename?  Load it

if nameOk % the file exists, read it.
    mIn=single(ReadEMFile(fullMergedImageName));
    ok=true;
else
    % try for reading the raw micrograph. We then subtract the median and
%     scale it to reflect fractional image intensity
    fullImageName=[mi.imagePath mi.imageFilenames{1}];
    if exist(fullImageName,'file') ...
            && isfield(mi,'imageNormScale') && mi.imageNormScale~=0;
        mIn=single(ReadEMFile(fullImageName));
        if ~isfield(mi,'imageMedian')
            mi.imageMedian=median(mIn(:));
        end;
        mIn=(mIn-mi.imageMedian)*mi.imageNormScale;
        ok=true;
    end;
end;
if ok
%     mIn=Crop(Downsample(mIn,mi.imageSize/2),mi.imageSize/2+4); %%%%%test code
    %     Determine what downsampling might already have occurred
    nIn=size(mIn); % we compare this with mi.imageSize
    M1=meGetImageScaling(mi.imageSize,nIn);

%     Determine our further downsampling
if numel(target)>1 % it's a target image size
        [M2,mOut]=GetImageScaling(mIn,target);
else % it's a maximum pixel size
        ds=floor(target/(M1(1,1)*mi.pixA)); % smallest integer downsampling that will work.
        nOut=NextNiceNumber(ceil(nIn/ds));
        [M2,mOut]=meGetImageScaling(mIn,nOut,ds);

end;
%  Combine the transformations.
M=M1*M2;
end;
