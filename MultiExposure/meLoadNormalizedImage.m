function [mOut,M,ok]=meLoadNormalizedImage(mi,targetSize,imgType)
% function [mOut,M,ok]=meLoadNormalizedImage(mi,targetSize,imgType)
%  -targetSize is a two-element vector, giving the dimensions of the
%  desired image, in pixels.
%  -M is the affine transform of mOut's possible downsampling and shift,
% Loads a merged image files, or the original micrograph
% assuming that we have normalization information in the mi structure.
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

maxScaleUp=1.1; % amount by which we'll allow a small image to be scaled up.
if numel(targetSize)<2 || targetSize(1)==0
    targetSize=mi.padImageSize;
end;
if nargin<3
    imgType='m'; % default: prefer small or compressed files.
end;
mergedImage=false;
rawImage=false;
M=zeros(3,3);
mOut=[];
ok=false;

if any(targetSize<mi.padImageSize) % Look for an existing small image, in the two possible places.
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
            m1=ReadEMFile(name);
            nIn=size(m1);
            if any(nIn*maxMag>=targetSize(1)) %% one dimension is big enough
                mergedImage=true;
            end;
        end;
        
    end;
end;
if ~mergedImage % try for a full-sized image
    name=[mi.procPath_sm mi.baseFilename imgType '.mrc'];
    if exist(name,'file')
        m1=ReadEMFile(name);
        nIn=size(m1);
        mergedImage=true;
    end;
end;

if ~mergedImage
    % try for reading the raw micrograph. We then subtract the median and
    %     scale it to reflect fractional image intensity
    fullImageName=[mi.imagePath mi.imageFilenames{1}];
    if exist(fullImageName,'file')
        m0=single(ReadEMFile(fullImageName));
        if ~all(size(m0)==mi.imageSize)
            error(['Micrograph size doesn''t match mi: ' num2str(size(m0)) ' vs ' num2str(mi.imageSize)]);
        end;
        nIn=mi.padImageSize;
        rawImage=true;
    end;
end;

ok=rawImage || mergedImage; % we got something
if ok
    dsMin=min(ceil(mi.padImageSize/targetSize(2))); % large dimension can't be bigger than maxTargetSize.
    dsMax=min(floor(mi.padImageSize/targetSize(1))); % large dimension can't be smaller than minTargetSize.
    ds=dsMin;
    for ds=dsMin:dsMax % no assignment if dsMin<1...
        if all(mod(nIn,ds)==0) % it's an integal divisor of the original micrograph
            break;
        end;
    end;
    if ds>0 % all set to downsample
        if rawImage % have to normalize and pad
            [mOut,M]=meMakeScaledMicrograph(m0,mi,ds);
        else
            M=meMakeMicrographScaleMatrix([],ds);
            nOut=mi.padImageSize/ds;
            ds1=nIn/ds;
            if any(mod(ds1,1)) % fractional downsampling
                mOut=DownsampleGeneral(mIn,nOut,1/ds1);
            else
                mOut=Downsample(m1,nOut);
            end;
        end;
    end;
end;