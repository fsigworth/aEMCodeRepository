function [effctz, mc, mts, coeffz, dctfs]=meCombineImages2(m,mi,ds,circMask,nZeros,mode)
% function [effctz, mc, mts, coeffz, dctfs]=meCombineImages2(m,mi,ds,circMask,nZeros,mode)
% Same as the old meCombineImages but requires that the mi structure is
% passed.  This allows movies to be supported.
% Now scales output images correctly, as identified by mi.mergeVersion=4
% Given a stack of original images m, the angstroms per pixel, the CTF
% parameters, and the array of composite transform matrices, combine the
% images into a single output image mc having the same noise level but
% optimal SSNR.  The effective CTF (including the ctf.ampFactors to model
% the apparent ctf) is returned in the 2D array effctz. The weights
% parameter provides an extra weighting of the individual images; by
% default it is all ones.
% - Individual aligned images are returned, if desired, in the stack mts.
% Those images are downsampled but are not filtered at all except by a circular
% mask if given.  To compute the overall merged image you would do
% mc = sum_i real( ifftn( fftn(mts(:,:,i)).* ifftshift(coeffz(:,:,i))))
% - The individual ctfs of the mts (with dose-dependent decay), and mi.weights
% not scaled by ctf.ampFactor are returned as dctfs.
% - The optional downsampling factor ds gives the
% "binning" of the output image; the default is 2.
% - The optional mode(1) parameter selects the final ctf form
%   1: normal, with the merging coefficients summing to 1 at each
%   frequency; this way the shot noise remains unchanged.  The merged image
%   is stored as *m.mrc.
%   2: simpleCTF, with merging done as above and then highpass filtered to
%   yield the same CTF as the low-defocus image but phase-flipped.  The net
%   result is a reduction of low-frequency noise.  The merged image is
%   stored as *msf.mrc
%   3: simpleCTF again, but with phases unflipped.  This is intended to mimick an
%   unprocessed micrograph, and is stored as *msu.mrc.
%   mode(2)=2 for old (pre-8/2016) merging algorithm; 3 for new one.
% set defaults
if nargin<3
    ds=2;
end;
if nargin<4
    circMask=0;
end;
if nargin<5
    nZeros=1;
end;
if nargin<6
    mode=1;
end;
if numel(mode)<2
    mode(2)=3;  % merging method 3, the latest
end;

CTFitPars=mi.ctf;
Tmats=mi.mergeMatrix;
doses=mi.doses;
if ~isfield(mi,'weights')
    mi.weights=ones(1,numel(doses));
end;
weights=mi.weights;
pixA=mi.pixA;

[nx, ny, nim]=size(m);
n0=[nx ny];
nd0=n0/ds;     % actual output image size;
if weights(1)>0 % we will use the first exposure, downsample everything.
    nd=nd0;     % the result had better be an integer!
else
    nd=2*nd0;  % we will use the later exposures only.  Oversample by 2 to
    % avoid interpolation artifacts on them
end;

for i=1:nim
    if ~isfield(CTFitPars,'ampFactor')
        CTFitPars(i).ampFactor=1;
    end;
end;

if nargout>2
    mts=zeros([nd0 nim],'single');  % Array to receive aligned images
end;

win=single(SquareWindow(n0,max(n0/512,4)));  % 8-point taper at 4k image

freqs=RadiusNorm(nd0)/(ds*pixA);  % output frequency samples

% Compute the mixing coefficients
if mode(2)==2
[coeffz, effctz, dctfs]=meComputeMergeCoeffs2(freqs, mi, nZeros, mode(1));
else
[coeffz, effctz, dctfs]=meComputeMergeCoeffs3(freqs, mi, nZeros, mode(1));
end;

% This code now inside meComputeMergeCoeffs2:
% switch mode
%     case 1 % use the coefficients as they are.  File *m.mrc
%     case 2 % reduce the lf coefficients.  Filename *msf.mrc
%         epsi=1e-6;
%         modFilter=abs(dctfs(:,:,1))./(effctz+epsi);
%         coeffz=coeffz.*repmat(modFilter,1,1,nim);
%         effctz=abs(dctfs(:,:,1));
%     case 3 % same, but restore alternating phases. Filename *msu.mrc
%         epsi=1e-6;
%         modFilter=dctfs(:,:,1)./(effctz+epsi);
%         coeffz=coeffz.*repmat(modFilter,1,1,nim);
%         effctz=dctfs(:,:,1);
%     otherwise
%         error('Unrecognized merging mode');
% end;

% Image sizes:
% n0 original micrograph size
% nd0 actual output merged image size = n0/ds
% nd intermediate output image size, 2 x nd0 if we upsample

if nargout >1  % Go ahead and compute the merged image
    mergeScale=1/ds^2;  % scale factor for downsampling
    %     Band-limit the first image and multiply by its coefficients
    if weights(1)>0
        if circMask
            fmskz=fuzzymask(nd,2,nd0*0.495,nd0*0.02);
        else
            fmskz=1;  % no mask
        end;
        fz0=Crop(fftshift(fftn(m(:,:,1))),nd).*fmskz*mergeScale;
        if nargout>2
            mts(:,:,1)=real(ifftn(ifftshift(fz0))); % Store the filtered image 1.
        end;
    else
        fz0=zeros(nd0);
    end;
    fz=fz0.*coeffz(:,:,1);
    % All the other images, which are being transformed,
    % will be reduced to final nyquist/2 bandwidth to avoid interp. artifacts
    fmskz=fuzzymask(nd,2,nd*0.245,nd*0.05); % prefilter for interpolation
    fmskz0=fuzzymask(nd0,2,nd0*0.45,nd0*0.05); % antialiasing for downsample
    %     T=Tmats(:,:,1);
    for i=2:nim
        if weights(i)>0
            fdz=Crop( fftshift( fftn(m(:,:,i).*win) ),nd).*fmskz*mergeScale;
            % Transform it and compute its ft
            md=AffineTransform(real( ifftn(ifftshift(fdz)) ),Tmats(:,:,i));
            if nd>nd0  % we have oversampled, downsample here to nd0
                md=real(ifftn(Cropo(fftn(md),nd0).*ifftshift(fmskz0)));
%                 Downsample(md,nd0,fmskz0); %% this would change scaling
            end;
            if nargout>2  % We return the transformed images
                mts(:,:,i)=md;
            end;
            fz0=coeffz(:,:,i).*fftshift(fftn(md));
            fz=fz+fz0;
        end;
    end;
    
    mc=real(ifftn(ifftshift(fz)));  % Return the combined image.
end;
