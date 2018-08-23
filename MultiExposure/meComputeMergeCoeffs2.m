function [coeffs, mergedCTF, dctfs, modFilter]=meComputeMergeCoeffs2( freqs, mi, nzeros, mode, noDamage )
% function [coeffs, mergedCTF, dctfs]=meComputeMergeCoeffs2( freqs, mi, nzeros, mode )
% Same as the old meComputeMergeCoeffs but it requires an mi structure as
% input, and supports damage compensation in movies.
%
% Given frequency values, mi.ctf an array of nims CTPars structures and a vector
% of nims doses, in e/A^2, compute the merging coefficients, which are the
% coefficients by which the fourier images are weighted before summing.
% Freqs (in A^-1) can be of any dimension. Suppose freqs is nx x ny.  Then
% the returned variable coeffs is nx x ny x nims.
% An example is freqs=RadiusNorm(n)/(ds*pixA);
% Image 1 is taken to be the lowest defocus.  The coeffs also effect phase
% flipping. The CTPars structure is the same as used for the
% ContrastTransfer() function.
% - doses is the dose in e/A^2.
% - The scalar or vector nzeros gives the
% number of zeros to preserve in a given image (image 1 always has the
% entire ctf preserved).  This prevents magnification of low-defocus images
% to fill zeros in high-defocus images. The default is 1.
% - mi.weights is an optional array of weights given the various images in
% computing coeffs.  To merge the first 2 of 3 images giving only .5
% amplitude to the second for example, you would set weights = [1 .5 0].
% - mergedCTF is the effective CTF computed according to the merging, and
% also includes the effect of radiation damage, the weights and
% ctf.ampFactors.  the latter account for the observation that the
% high-defocus images actually contain more signal than expected from the
% theoretical ctfs.  So we don't use the ampFactors in computing the coeffs
% but we do put them into the mergedCTF.
% - dctfs is an array of the 'raw' ctfs of each image.  They include radiation
% damage and are truncated after nzeros zeros.  They are not multiplied by
% weights or ctf.ampFactors.  These are the ctfs used to compute the
% coefficients.  Multiply them by doses(i) * ctf(i).ampFactor *
% coeffs(:,:,i), and sum them to get the overall effective ctf the merged
% image, = mergedCTF which is returned.
%
% - The optional mode parameter selects the final ctf form
%   1: normal, with the merging coefficients summing to 1 at each
%   frequency; this way the shot noise remains unchanged.  The merged image
%   is stored as *m.mrc.
%   2: simpleCTF, with merging done as above and then highpass filtered to
%   yield the same CTF as the low-defocus image but phase-flipped.  The net
%   result is a reduction of low-frequency noise.  The merged image is
%   stored as *msf.mrc
%   3: simpleCTF again, but with phases unflipped.  This is intended to mimick an
%   unprocessed micrograph, and is stored as *msu.mrc.
%   4: normal, but with no damage compensation
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
epsi=1e-12;  % smallest noise variance to preserve.
% underflow of single.

if nargin<5
    noDamage=0;
end;
if nargin<4
    mode=1;
end;

if nargin<3
    nzeros=1;
end;

% Pick up parameters from the mi structure
CTPars=mi.ctf;
doses=mi.doses;
weights=mi.weights;
if ~isfield(mi.ctf(1),'ampFactor')
    for i=1:numel(CTPars)
        mi.ctf(i).ampFactor=1;
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

% Let the first exposure be the standard
dose1=doses(1);

%%
% Compute the decay constant due to radiation, as a function of frequency
% --a fit to the data in the Rubenstein paper.
%%%%%% n0=2*0.16./(abs(f)+.004)+5;  % twice the critical dose at 200 kV
% modified to better match vesicle decay.
% n0=.32./(abs(f.^2*10)+.002)+5;
if isfield(mi,'damageModelCode') && numel(mi.damageModelCode)>1
    n0=eval(mi.damageModelCode);
else
    if ~noDamage
        warning('No damageModelCode in the mi structure');
    end;
    mi.damageModelCode = '.184*f.^-1.665+2.1;';
    mi.weights=[1 1];
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
% There are two possibilities: one continuous movie with sequenctial
% frame sets, or two or more movies.  We identify multiple movies from the
% class and number of elements of mi.movieFilename
pastAccumDose=0;  % remember the dose from one movie to the next
for i=find(doses(:)'>0)
    if noDamage
        signal=1;
    else
        % compute the normalized signal amplitude after decay (=1
        %     given no decay)
        if isfield(mi,'frameDose') && numel(mi.frameDose)>0 ...
                && all(mi.frameDose>0)  % we have a movie: do a calculation like
            %         in k2DamageCompCTF.
            startFrame=mi.frameSets(i,1);
            endFrame=mi.frameSets(i,2);
            nFrames=endFrame-startFrame+1;
            if i>1 && mi.frameSets(i,1)>mi.frameSets(i-1,2) % a continuous movie
                pastAccumDose=0;
            end;
            %         mi.damageModelCode='2*(.245*abs(f).^-1.665+2.81)';
            %         if ~isfield(mi,'damageModelCode')
            %             mi.damageModelCode='0.32./(abs(f.^2*10)+.002)+5;';
            %             warning('Used the default damage model code');
            %         end;
            n0=min(1e5,eval(mi.damageModelCode));  % a function of f
            k0=1./n0;  % decay rate
            % Make the raw transfer functions as sqrt(SNR)
            if numel(mi.frameDose)<2
                d=mi.frameDose*ones(nFrames,1);
            else
                d=mi.frameDose(startFrame:endFrame);
            end;
            % ? are we perhaps skipping dose before the startFrame??
            accumDose=pastAccumDose+[0 cumsum(d(:))'];
            G=zeros(nel,1,'single');
            %         There is a more efficient way to do this, using geometric
            %         series, but here is brute force...
            for j=1:nFrames
                Q=1-exp(-d(j)*k0);
                h=(n0/d(j)).*exp(-accumDose(j)*k0).*Q;
                G=G+h.^2;
            end;
            % normalize to give white shot noise
            signal=sqrt(G/nFrames);
            pastAccumDose=accumDose(end);
        else
            signal=(exp(-cumdose(i)./n0).*n0.*(1-exp(-doses(i)./n0)))/doses(i);
        end;
    end;
    [ctf,chi]=ContrastTransfer(freqs,CTPars(i));
    fmask=abs(chi(:))<nzeros(i);  % frequency mask beyond zeros.
    %     ctf including decay, mask and amp factor
    dctfs(:,i)=ctf(:).*signal.*fmask;  % no ampFactor or weight is applied.
    %     coeffs reflect the weights.
    coeffs(:,i)=dctfs(:,i)*weights(i)*doses(i)/dose1;  % sum of dose-weighted ctfs
end;
coeffs(isnan(coeffs))=0;
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
coeffs(isnan(coeffs))=0;

dctfs=reshape(dctfs,sizx);
mergedCTF=reshape(mergedCTF,siz);

% Modify the coeffs and mergedCTF according to the merging mode
switch mode
    case 1 % use the coefficients as they are.
        modFilter=1;
    case 2 % reduce the lf coefficients.
        epsi=1e-6;
        modFilter=abs(dctfs(:,:,1))./(mergedCTF+epsi);
        coeffs=coeffs.*repmat(modFilter,1,1,nim);
        mergedCTF=abs(dctfs(:,:,1));
    case 3 % same, but restore alternating phases.  No inversion of lowest frequencies.
        epsi=1e-6;
        modFilter=-dctfs(:,:,1)./(mergedCTF+epsi);
        coeffs=coeffs.*repmat(modFilter,1,1,nim);
        mergedCTF=-dctfs(:,:,1);
    case 4 % use the coefficients as they are
    otherwise
        error('Unrecognized merging mode');
end;

mergedCTF(isnan(mergedCTF))=0;



% freqs=reshape(freqs,siz);

% % more test code...
% subplot(2,1,1);
% semilogx(freqs,[abs(signals.*ctfs)/10 mergedCTF]);
%
% subplot(2,1,2);
% semilogx(freqs,abs(coeffs));
%
% % imacs(mergedCTF);
