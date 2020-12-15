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

maxScaleUp=1.2; % amount by which we'll allow a small image to be scaled up.
if numel(targetSize)<2 || targetSize(1)==0
    targetSize=mi.padImageSize;
end;
if nargin<3
    imgType='m'; % default: prefer small or compressed files.
end;
ok=false;
M=zeros(3,3);
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
        nIn=size(mIn);
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
        nIn=mi.padImageSize;
        mIn=Crop(m0,mi.padImageSize);
        ok=true;
    end;
end;

if ok
    ds1=max(mi.padImageSize./nIn);
    M1=meMakeMicrographScaleMatrix(mi,ds1);  % matrix to produce the image
%     we've got so far.
    ds2=max(nIn./targetSize); % we pick the dimension best matched by the target.
    mOut=DownsampleGeneral(mIn,targetSize,1/ds2);
    shift=floor(nIn./ds2-targetSize)/2; % The shift from the crop operation
    M2=[ds2 0 -shift(1); 0 ds2 -shift(2); 0 0 1];
    M=M1*M2; % Composite matrix mapping output to original image
end;