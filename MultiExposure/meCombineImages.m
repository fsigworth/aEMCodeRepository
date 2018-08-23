function [effctz, mc, mts, coeffz, dctfs]=meCombineImages(m,pixA,CTFitPars,Tmats,doses,ds,circMask,weights,nZeros)
% function [effctz mc mts coeffz dctfs]=meCombineImages(m,pixA,CTFitPars,Tmats,doses,ds,circMask,weights,nZeros)
% function [effctz mc mts coeffz dctfs]=meCombineImages(m,mi,ds,circMask,nZeros)
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
% Changed to use meComputeMergeCoeffs2, which takes mi as the main
% argument.

% Pick up either 9 or 5 arguments, and set defaults
if isa(pixA,'struct')
    mi=pixA;
    if nargin>=3
        ds=CTFitPars;
    else
        ds=2;
    end;
    if nargin>=4
        circMask=Tmats;
    else
        circMask=0;
    end;
    if nargin>=5
        nZeros=doses;
    else
        nZeros=1;
    end;
    CTFitPars=mi.ctf;
    Tmats=mi.mergeMatrix;
    doses=mi.doses;
    if ~isfield(mi,'weights')
        mi.weights=ones(1,numel(doses));
    end;
    weights=mi.weights;
    pixA=mi.pixA;
else
    
    mi.pixA=pixA;
    mi.ctf=CTFitPars;
    mi.doses=doses;
    mi.mergeMatrix=Tmats;
    
    if nargin<6
        ds=2;  % default is to downsample by 2.
    end;
    if nargin<7
        circMask=0;  % Make a circular mask in Fourier space.
    end;
    if nargin<8
        weights=ones(1,numel(doses));
    end;
    if nargin<9
        nZeros=1;
    end;
    mi.weights=weights;
end;


[nx ny nim]=size(m);
n0=[nx ny];
nd0=n0/ds;     % actual output image size;
if weights(1)>0 % we will use the first exposure, downsample everything.
    nd=nd0;     % the result had better be an integer!
else
    nd=2*nd0;  % we will use the later exposures only.  Oversample by 2 to
                % avoid interpolation artifacts on them
end;

for i=1:3
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
% [coeffz, effctz, dctfs]=meComputeMergeCoeffs(freqs, CTFitPars, doses, nZeros, weights);
[coeffz, effctz, dctfs]=meComputeMergeCoeffs2(freqs, mi, nZeros);

if nargout >1  % Go ahead and compute the merged image
    %     Band-limit the first image and multiply by its coefficients
    if weights(1)>0
        if circMask
            fmskz=fuzzymask(nd,2,nd*0.495,nd*0.02);
        else
            fmskz=1;  % no mask
        end;
        fz0=Crop(fftshift(fftn(m(:,:,1))),nd).*fmskz;
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
        fdz=Crop( fftshift( fftn(m(:,:,i).*win) ),nd).*fmskz;
        % Transform it and compute its ft
        md=AffineTransform(real( ifftn(ifftshift(fdz)) ),Tmats(:,:,i));
        if nd>nd0  % we have oversampled, downsample here to nd0
            md=Downsample(md,nd0,fmskz0);
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
