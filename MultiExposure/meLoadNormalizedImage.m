function [mOut,M,ok,rawImg]=meLoadNormalizedImage(mi,targetSize,imgType)
% function [mOut,M,ok]=meLoadNormalizedImage(mi,targetSize,imgType)
%  -targetSize is a two-element vector, giving the dimensions of the
%  desired image, in pixels.
%  -M is the affine transform of mOut's possible downsampling and shift,
% Loads a merged image files, or the original micrograph
% assuming that we have normalization information in the mi structure.
% rawImg=true if we loaded the raw micrograph, false if we loaded some sort
% of merged (i.e. padded) image, even though we convert it to a padded
% image anyway.
% We assume that the size of any 'small' merged image we read will be a submultiple
% of mi.padImageSize, i.e. has an integer downsampling factor.
% 
% imgType is a string such as 'm' or 'mv', used to construct filenames like
% xxxxmvs.mrc or xxxxmv.mrc
% Return an image mOut that is cropped and downsampled to be
% targetSize and have an integer downsampling factor from
% the original micrograph. So if the original image is rectangular but
% targetSize is square, a square output image is provided, but with
% appropriate offsets put into M.
% The returned matrix M describes the coordinate transformation with
% scaling and possible padding.
%   M(1,1): downsampling factor from original micrograph to mOut
%   -M(1:2,3): shift added to original micrograph before downsampling. Thus to
%   determine a pixel coordinate in the the original micrograph, it is
%   computed so:
% origMicrographXY=M*[xOut;yOut;1];
% [xOut;yOut;1]=M\[origX origY 1];
%  assuming zero-based coordinates in each case.

maxScaleUp=1.2; % amount by which we'll allow a small image to be scaled up.
if nargin<2
    targetSize=mi.padImageSize;
elseif numel(targetSize)<2
    targetSize=targetSize*[1 1];
end;
if nargin<3
    imgType='m'; % default: prefer small or compressed files.
end;
ok=false;
rawImg=false;
M=eye(3); % default is an identity transformation
origOffset=-floor((mi.padImageSize-mi.imageSize)/2);

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
    name=[mi.procPath mi.baseFilename imgType '.mrc'];
    if exist(name,'file')
        mIn=ReadEMFile(name);
        ok=true;
    end;
end;

if ~ok && ~any(imgType=='v')
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
    ds=ceil(max(mi.padImageSize./targetSize)); % desured overall downsampling factor.
    M(1:2,3)=(origOffset-floor((ds*targetSize-mi.padImageSize)/2))';
    M(1,1)=ds;
    M(2,2)=ds;
    ds1=min(mi.padImageSize./size(mIn)); % initial downsampling
    mOut=DownsampleGeneral(mIn,targetSize,ds1/ds); % additional downsampling
end;
