function [mOut,M,ok,rawImg]=meLoadNormalizedImage(mi,targetSize,imgType)
% function [mOut,M,ok]=meLoadNormalizedImage(mi,targetSize,imgType)
%  -targetSize is a two-element vector, giving the dimensions of the
%  desired image, in pixels.
%  -M is the affine transform of mOut's possible downsampling and shift,
% Loads a merged image files, or the original micrograph
% assuming that we have normalization information in the mi structure.
% rawImg=true if we loaded the raw micrograph, false if we loaded some sort
% of merged (i.e. padded) image.
% 
% imgType is a string such as 'm' or 'mv', used to construct filenames like
% xxxxmvs.mrc or xxxxmv.mrc
% Return an image mOut that is cropped and downsampled to be at most
% targetSize and have an integer downsampling factor from
% the original micrograph.
% The returned matrix M describes the coordinate transformation with
% scaling and possible padding.
%   M(1,1): downsampling factor from original micrograph to mOut
%   -M(1:2,3): shift added to original micrograph before downsampling. Thus to
%   determine a pixel coordinate in the the original micrograph, it is
%   computed so:
% origMicrographXY=M*[xOut;yOut;1];
%  assuming zero-based coordinates.

maxScaleUp=1.2; % amount by which we'll allow a small image to be scaled up.
if numel(targetSize)<2 || targetSize(1)==0
    targetSize=mi.padImageSize;
end;
if nargin<3
    imgType='m'; % default: prefer small or compressed files.
end;
ok=false;
rawImg=false;
M=zeros(3,3);

M1=eye(3); % identity transformation
M1(1:2,3)=-floor((mi.padImageSize-mi.imageSize)/2);

mOut=[];

if any(targetSize<mi.padImageSize) % Look for an existing small image *s.mrc,
    %     in the two possible places.
    paths=cell(1);
    if isfield(mi,'procPath_sm') % First look here, if the directory exists
        paths={mi.procPath_sm};
    end;
    paths=[paths {mi.procPath}]; % Old pattern, where *ms.mrc or *mvs.mrc
    %     files were stored.
    
    for iPath=1:numel(paths)
        % Load a small image
        name=[paths{iPath} mi.baseFilename imgType 's.mrc'];
        if exist(name,'file')
            mIn=ReadEMFile(name);
            nIn=size(mIn);
            if all(nIn*maxScaleUp>=targetSize(1)) %% it's big enough
                ok=true;
                break;
            end;
        end;
        
    end;
end;
if ~ok % try for a full-sized image in the procPath directory
    name=[mi.procPath_sm mi.baseFilename imgType '.mrc'];
    if exist(name,'file')
        mIn=ReadEMFile(name);
        ok=true;
    end;
end;

if ~ok
    % try for reading the raw micrograph. We then subtract the median and
    %     scale it to reflect fractional image intensity, and pad it.
    fullImageName=[mi.imagePath mi.imageFilenames{1}];
    if exist(fullImageName,'file')
        m0=single(ReadEMFile(fullImageName));
        if ~all(size(m0)==mi.imageSize)
            error(['Micrograph size doesn''t match mi: ' num2str(size(m0)) ' vs ' num2str(mi.imageSize)]);
        end;
        mIn=Crop((m0-mi.imageMedian)*mi.imageNormScale,mi.padImageSize);
        rawImg=true;
        ok=true;
    end;
end;


if ok
%     Fill in the scaling of the image we read in.
    ds=mi.padImageSize./size(mIn); % downsampling factor of image we're given.
    M1(1,1)=1/ds(1);
    M1(2,2)=1/ds(2);
    
    [mOut,M]=meDownsampleImage(mIn,M1,targetSize);
end;
