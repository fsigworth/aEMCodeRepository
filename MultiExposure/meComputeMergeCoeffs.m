function [coeffs, mergedCTF, dctfs]=meComputeMergeCoeffs( freqs, CTPars, doses, nzeros, weights)
% [coeffs mergedCTF dctfs]=meComputeMergeCoeffs(freqs,CTPars,doses,nzeros,weights)
% or: [coeffs mergedCTF dctfs]=meComputeMergeCoeffs( freqs, mi, nzeros)
% obsolete version, use meComputeMergeCoeffs2.
% Given frequency values, an array of nims CTPars structures and a vector
% of nims doses, in e/A^2, compute the merging coefficients, which are the
% coefficients by which the fourier images are weighted before summing.
% Freqs (in A^-1) can be of any dimension. Suppose freqs is nx x ny.  Then
% the returned variable coeffs is nx x ny x nims.
% Image 1 is taken to be the lowest defocus.  The coeffs also effect phase
% flipping. The CTPars structure is the same as used for the
% ContrastTransfer() function.
% - doses is the dose in e/A^2. The scalar or vector nzeros gives the
% number of zeros to preserve in a given image (image 1 always has the
% entire ctf preserved).  This prevents magnification of low-defocus images
% to fill zeros in high-defocus images. The default is 1.
% - weights is an optional array of weights given the various images in
% computing coeffs.  To merge the first 2 of 3 images giving only .5
% amplitude to the second for example, you would set weights = [1 .5 0].
% - mergedCTF is the effective CTF computed according to the merging, and
% also includes the effect of radiation damage, the weights and
% ctf.ampFactors.  the latter account for the observation that the
% high-defocus images actually contain more signal than expected from the
% theoretical ctfs.  So we don't use the ampFactors in computing the coeffs
% but we do put them into the mergedCTF.
% - dctfs is an array of the 'raw' ctfs of each image, including radiation
% damage and truncated after nzeros zeros.  They are not multiplied by
% weights or ctf.ampFactors.  These are the ctfs used to compute the
% coefficients.  Multiply them by doses(i) * ctf(i).ampFactor *
% coeffs(:,:,i), and sum them to get the overall effective ctf the merged
% image, = mergedCTF which is returned.
%
% Example 1, to graph the effect of merging,
% f=0:.001:.1;
% [coeffs mergedCTF]=ComputeMergeCoeffs(f,CTPars,[10 20]);
% plot(f, mergedCTF);
%
% Example 2, to do merging of images m1 and m2, already aligned:
% f1=fftn(m1);
% f2=fftn(m2);
% freq=fftshift(Radius(n)/(n*pixA));  % image size is n x n, and the
%                                     % pixel size is pixA angstroms.
% [coeffs ceff]=ComputeMergeCoeffs( freq, CTPars, [10 20]);
% fmerge=f1.*coeffs(:,:,1)+f2.*coeffs(:,:,2);
% mmerge=real(ifftn(fmerge));
%
% memory use for a 4k x 4k frequency array is about (6+3*nim)*68MB.
% fixed modeling of decay (using dctfs in computing merged CTF) fs Apr 2013

% % test code:
% nargin=5;
% freqs=(0:.0001:.3)';
% % freqs=Radius(4096)/1e4;
% % disp('start');
% %
% P.lambda=.025;
% P.Cs=2;
% P.alpha=.07;
%
% defvals=[.6 3 10];
% doses=[10 10 20];
% nzeros=1;
% weights=[1 1 1];
% for i=1:3
%     P.defocus=defvals(i);
%     P.B=100*defvals(i);
%     CTPars(i)=P;
% end;

% -----------
epsi=1e-20;  % smallest noise variance to preserve.
kV=200;  % assumed energy.

% underflow of single.

% Handle the two options for arguments.
if isfield(CTPars,'ctf') % This is actually an mi structure
    mi=CTPars;
    CTPars=mi.ctf;
    nzeros=doses;
    doses=mi.doses;
    weights=mi.weights;
    if ~isfield(mi.ctf(1),'ampFactor')
        for i=1:numel(CTPars)
            mi.ctf(i).ampFactor=1;
        end;
    end;
    if nargin<3
        nzeros=1;
    end;
else
    if nargin<4
        nzeros=1;
    end;
    if nargin<5
        weights=ones(1,numel(CTPars));
    end;
end;


nim=numel(CTPars);

if numel(nzeros)<nim
    nzeros=ones(nim,1)*nzeros;
end;

weights=weights(:);

siz=size(freqs);

% create the dimension of the output variables
sizx=siz;
ndims=numel(sizx);
if sizx(ndims)==1  % remove a trailing singleton dimension
    sizx=siz(1:ndims-1);
end;
sizx=[sizx nim];

nel=numel(freqs);  % 1d variable index

f=abs(single(freqs(:)));  % take the absolute value!
doses=single(doses(:));
cumdose=cumsum([0; doses]);

% find the first exposure with a nonzero dose, and let that be the standard
q=find(doses,1);
dose1=doses(q);

% Compute the decay constant due to radiation, as a function of frequency
% --a fit to the data in the Rubenstein paper.
%%%%%% n0=2*0.16./(abs(f)+.004)+5;  % twice the critical dose at 200 kV
% modified to better match vesicle decay.
    if kV>200
        n0=.245*f.^-1.665+2.81; % Grant&Grigorieff 300 kV'
    else
        n0=.184*f.^-1.665+2.1; % Grant&Grigorieff 200 kV'
    end;


iFirst=find(weights,1); % mark the first nonzero weight
nzeros(iFirst)=inf;  % by default, no masking of 1st exposure.

% Get the coefficients coeff = ctf*signal
% Coefficients should be proportional to signal amp / rms noise
% where signal is the relative signal after radiation damage.
dctfs=single(zeros(nel,nim));
coeffs=single(zeros(nel,nim));

if ~isfield(CTPars(nim),'ampFactor') || numel(CTPars(nim).ampFactor)<1
    for i=1:nim
        CTPars(i).ampFactor=1;
    end;
end;
for i=find(doses(:)'>0)
    % compute the normalized signal amplitude after decay (=1
    %     given no decay)
    if isfield(mi,'frameDose') && numel(mi.frameDose)==1 && mi.frameDose>0  % we have a movie: do a calculation like
%         in k2DamageCompCTF.
        startFrame=mi.frameSets(segIndex,1);
        nFrames=mi.frameSets(segIndex,2)-startFrame+1;
        n0=eval(mi.damageModelCode);  % a function of f
        k0=1./n0;  % decay rate
        % Make the raw transfer functions as sqrt(SNR)
        d=mi.frameDose;
        G=zeros(nel,1,'single');
        Q=1-exp(-d*k0);
%         There is a more efficient way to do this, using geometric
%         series, but brute force...
        for j=1:nFrames
            dStart=d*(startFrame+j-2);
            h=(n0/d).*exp(-dStart*k0).*Q;
            G=G+h.^2;
        end;
        % normalize to give white shot noise
        signal=sqrt(G/nFrames);
    else
        signal=(exp(-cumdose(i)./n0).*n0.*(1-exp(-doses(i)./n0)))/doses(i);
    end;
    [ctf,chi]=ContrastTransfer(freqs,CTPars(i));
    fmask=abs(chi(:))<nzeros(i);  % frequency mask beyond zeros.
    %     ctf including decay, mask and amp factor
    dctfs(:,i)=ctf(:).*signal.*fmask;  % no ampFactor or weight is applied.
    %     coeffs reflect the weights.
    coeffs(:,i)=dctfs(:,i)*weights(i)*doses(i)/dose1;  % sum of dose-weighted ctfs
end;

% compute the frequency-dependent denominator
D=single(zeros(nel,1));
for i=1:nim
    D=D+coeffs(:,i).^2;  % sum the noise variance
end;
D(D<epsi)=1;
D=sqrt(D);  % noise standard deviation given non-normalized coefficients
%%
mergedCTF=single(zeros(nel,1));
for i=1:nim
    coeffs(:,i)=coeffs(:,i)./D;  % normalize to signal at dose=1.
    %     Note that in computing the MergedCTF we *do* use the ampFactors.
    mergedCTF=mergedCTF+coeffs(:,i).*dctfs(:,i)*CTPars(i).ampFactor*doses(i)/dose1;
end;

coeffs=reshape(coeffs,sizx);
dctfs=reshape(dctfs,sizx);
mergedCTF=reshape(mergedCTF,siz);

% freqs=reshape(freqs,siz);

% % more test code...
% subplot(2,1,1);
% semilogx(freqs,[abs(signals.*ctfs)/10 mergedCTF]);
%
% subplot(2,1,2);
% semilogx(freqs,abs(coeffs));
%
% % imacs(mergedCTF);
