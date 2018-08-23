function roi=reEMStep26(imgs,refs,si,ri,moi,altImgs,logs)
% function roi=reEMStep26(imgs,refs,si,ri,moi,altImgs)
% Like reEMStep24, but reduces memory requirement by operating on one
% micrograph's images at a time.
% Like reEMStep23, but puts all the output variables into roi, and uses
% ri.flags fields to determine which fields are computed.  Also
% moi.activeFlags and roi.pTrans are both nt^2 x nim, instead of nt x nt x
% nim as before.
% imgs is either n x n x nImgs or n x n x nAlphas*nIsos x nImgs, pre-rotated.
% if present, altImgs must be the same size.
% refs is a corresponding stack n x n x nRefs x nVols in size.
% si is the stack info structure, from which we use ctfs, miIndex,
%   rVesicle, yClick, mbnOffset.
% ri is the run info which contains
%   nTrans, alphas, isos, radius, softMask, angleStep andgleMin, angleN,
%   useAltImgs
% moi is the model parameters, including sigmaN,sigmaC,sigmaG,pVols,
%   and imgAmps (nImgs x 1) amd also pRefs, which is either a scalar
%   or nImgs x nRefs x nVols in size.
% moi.activeTrans is nt^2 x nImgs booleans indicating which translations will be
%   pursued.
% altImgs is an optional stack of images (e.g. membrane not subtracted) the
% same size as imgs but used for the M-step only.  Thus the latent
% probabilities are computed from imgs, but altImgs are used to for the
% class means.
%
% The run output is provided in the structure roi, see the end of this
%   file for its fields.
% imgBestMatch is the transformed, marginalized over alpha and translation,
% image corresponding to the best reference given by roi.imgClass.
%
% This is a fast version that uses brute-force translations.  On my 4-core
% machine It runs faster with 'for' rather than 'parfor' loops presumably
% because of Matlab's multithreaded FFT and linear algebra functions.
%
% ri.radius (upper bound on hard mask radius), ri.softMask
% (all are assumed to have the same size, nw = ri.n, as imgs and refs)
useAltImgs=(nargin>5) && (numel(altImgs)>0);
makeRawProbs=isfield(ri.flags,'makeRawProbs') && ri.flags.makeRawProbs;

sz=size(imgs);
nw=sz(1);     % Working image size
preRotated=numel(sz)>3;
nImgs=sz(end);

inBoundImgs=true(1,nImgs);  % booleans keeping track of images that were included in means

nt=ri.nTrans;    % Size of translations calculation, must be odd
nHalfTrans=floor((nt-1)/2);
nt2=nt^2;               % total trans. pixels to calculate
nT=nt2;

nMi=size(si.ctfs,3);    % number of micrographs

nRefs=size(refs,3);
nVols=ri.nVols;

alphasI=ri.alphas(:);  % list of all alphas including isos
nAlphas=size(ri.alphas,1);
nGammas=ri.angleN(3);
nBetas=ri.angleN(2);
nIsos=size(ri.alphas,2);

% Define bad-particle values for alphaI and translations
badAI=false(size(ri.alphas)); % extreme alpha values are bad
badAI(1,:)=true;
badAI(end,:)=true;
badTy=[1 nt];
badTx=[1 nt];
%badTx=[1:1+floor(nHalfTrans*(1-ri.maxTxFraction)) 1+ceil(nHalfTrans*(1+ri.maxTxFraction)):nt];

nRV=nRefs*nVols;
nAI=nAlphas*nIsos;

% Priors on volumes and references
logPVols=max(log(moi.pVols),-12);  % probs are > 6e-6.
logPRVs=max(log(moi.pRefs(:)),-12); % pRefs is nRefs x nVols
if numel(logPRVs)==1
    logPRVs=logPRVs*ones(nRV,1);
end;

% Create the hard mask for selecting pixels
% pick the radius for the hard mask to be maximum to allow translations
hmRadius=floor(min(ri.radius,floor(nw-1)/2-nHalfTrans));


% Mask for centered particles
msk=fuzzymask(nw,2,hmRadius,0); % r should be integer, 1 more than actual radius
invSoftMask=1-ri.softMask;   % annular mask
invSoftMask=invSoftMask/sum(invSoftMask(:)); % sums to 1
nPixEff=ri.softMask(:)'*ri.softMask(:); % softMask has pixels approaching 1

pix=msk(:);  % nw^2 x 1: picks out the pixels for centered mask
nPix=sum(pix);
softMaskPix=ri.softMask(pix);

% For particles with all shifts
transPix=false(nw^2,nt^2);
offset=floor(nw/2+1)-nHalfTrans-1;

% mask is translated up and right with increasing ix,iy, so
% img(transPix) is translated left and down.  So if particle is in upper
% right, positive ix,iy brings it to the center.
for ix=1:nt
    for iy=1:nt
        msk=fuzzymask(nw,2,hmRadius,0,[ix iy]+offset);
        transPix(:,ix+(iy-1)*nt)=msk(:);
    end;
end;
% Check for consistency
if ~all(sum(transPix,1)==sum(pix))
    error('Hard mask size varies with translations');
end;

% Allocate the output variables
pAlphas=zeros(nAlphas*nIsos,nImgs,'single');
pRefs=zeros(nRV,nImgs,'single');

% Allocate the other output variables
logPX=zeros(1,nImgs);  % Likelihood.  Note that unassigned values are zero!
imgAmps=zeros(1,nImgs,'single');
eVar=zeros(1,nImgs,'single');  % noise variance
pT_X0s=zeros(nt2,nImgs,'single');  % prob of (absolute) translations
bmClass=zeros(1,nImgs,'single');  % best-matching ref/volume
bmTAI=zeros(3,nImgs,'single'); % best matching trans/alpha/iso
pixBestMatch=zeros(nPix,nImgs,'single'); % best transformed image

%     Accumulate the classMeans and classNorms using the ctfs
fclassMeans=zeros(nw,nw,nRV,'single');  % will be complex
classNorms=zeros(nw*nw,nRV,'single');

%********
if makeRawProbs
    rawLogProbs.priors=zeros(nT,nAI,nRV,nImgs,'single');
    rawLogProbs.pXs=zeros(nT,nAI,nRV,nImgs,'single');
end;
%*********

% Sort all the particles according to micrograph
miPartIndices=cell(nMi,1); % cell array of row vectors of image indices
for i=1:nImgs
    miPartIndices{si.miIndex(i)}(1,end+1)=i;
end;


% Enumerate all the alphas including inside-out flips
% alphasI=zeros(nAlphas,nIsos,'single');
% for isi=1:nIsos
%     io=ri.isos(isi);
%     alphasI(:,isi)=alphas(:)+180*io;
% end;
% alphasI=alphasI(:)';

%     create the rotation matrix for translation
% rotAlphaMat=reMakeTransAlphaMatrix(alphasI,nt);

% -------------------Loop over micrographs------------------------
imageCount=0;
totalCount=0;
tic
for iMi=1:nMi
    % Variables accumulated for each micrograph
    pixMeans=zeros(nPix,nRV);
    rawNorms=zeros(nRV,1);
    ctRefs=zeros(nPix,nRefs,nVols);
    R2=zeros(nRefs,nVols);  % 2nd moment of refs
    
    ctf=si.ctfs(:,:,iMi);
    ctfo=ifftshift(ctf);
    %    create the references, filtered by ctfs.
    %     They are multiplied by softMask and sampled into pixels
    for iVol=1:nVols
        for iRef=1:nRefs
            fref=fftn(refs(:,:,iRef,iVol));
            ctRef0=real(ifftn(fref.*ctfo));
            dc=invSoftMask(:)'*ctRef0(:);
            ctRef=ri.softMask.*(ctRef0-dc);
            ctRefPix=ctRef(pix);
            R2(iRef,iVol)=ctRefPix'*ctRefPix;
            ctRefs(:,iRef,iVol)=ctRefPix;
        end;
    end;
ctRefs=reshape(ctRefs,nPix,nRV);
R2=reshape(R2,nRV,1);
    
    % ---------Loop over images in a micrograph---------
    for iImg=miPartIndices{iMi} % all the particles from this micrograph
        % write something every 5 min or 1/4 of set to show we're alive
        imageCount=imageCount+1;
        totalCount=totalCount+1;
        ti=toc;
        if ti>300 || imageCount>nImgs/4
            mdisp(logs,['  ' datestr(now) ' image ' num2str(totalCount)]);
            imageCount=0;
            tic;
        end;
        
        imgAmp=moi.imgAmps(iImg);
        tc=moi.activeTrans(:,iImg);
        aic=moi.activeAlphas(:,:,iImg);
        if imgAmp>0 && any(tc(:)) && any(aic(:))  % This image is active
            % Make the pointers to sparse translations
            tc=tc(:);  % active flags: booleans
            tp=find(tc);  % pointers of the active -> all translations
            ntc=numel(tp);  % number of active translations
%             Make pointers to sparse alphas and isos
            aip=find(aic(:));
            nai=numel(aip);
            ntai=ntc*nai;  % lower case means active values only.
%               Sparse references
            arv=moi.activeRVs(:,:,:,iImg); % boolean nRefs x nVols
            arv=reshape(arv,nRefs,nVols);
            narv=sum(arv(:));
            rvp=find(arv(:));
            aVols=floor((rvp-1)/nRefs)+1; % vol index for each arv

%             Select the references and ref powers
            ctaRefs=ctRefs(:,arv(:));
            Ra2=R2(arv(:));
            
            vesR=si.rVesicle(iImg);
            y0=si.yClick(iImg);
            
            %         ------set up the images------
            %   Images are rotated ccw.  That means if particle is rotated to
            %   right (cw rotated) then a positive alpha is the one that brings the particle into
            %   register.
            if preRotated
                unrotImgs=imgs(:,:,aip,iImg);
                %           Get the power from an image with angle near 0
                avals=mod(alphasI(aip),90);
                avals(avals>45)=avals(avals>45)-90;
                [~,iAlpha0]=min(abs(avals));
                img=unrotImgs(:,:,iAlpha0); % pick this un-rotated one for power calc.
            else
                img=imgs(:,:,iImg);
                unrotImgs=rsRotateImage(img,-alphasI(aip));
            end;
%         Create the matrix of shifted and rotated images
            transImgs=zeros(nPix,ntc,nai);  % 1D images nTrans*nAlpha
            %         Make the translated images, multipled by soft mask
            for i=1:nai  % loop over alphas
                uImg=unrotImgs(:,:,i);
                uImg=uImg(:);
                for j=1:ntc  % loop over translations
                    transImgs(:,j,i)=uImg(transPix(:,tp(j))).*softMaskPix;
                end;
            end;
            transImgs=reshape(transImgs,nPix,ntai);  % trans change most rapidly.
            
            %         ------set up for altImgs------
            if useAltImgs  % in case we're using alternative images for recon
                if preRotated
                    unrotAltImgs=altImgs(:,:,aip,iImg);
                else
                    unrotAltImgs=rsRotateImage(altImgs(:,:,iImg),-alphasI(aip));
                end;
                transAltImgs=zeros(nPix,ntc,nai);
                for i=1:nai  % loop over alphas
                    uAltImg=unrotAltImgs(:,:,i);
                    uAltImg=uAltImg(:);
                    for j=1:ntc  % loop over translations
                        transAltImgs(:,j,i)=uAltImg(transPix(:,tp(j))).*softMaskPix;
                    end;
                end;
                transAltImgs=reshape(transAltImgs,nPix,ntai);  % trans change most rapidly.
            end;
            
            % -------- Here's the big matrix multiplication --------
            ccIR=transImgs'*ctaRefs;  % cross-correlation, ntAI x nRV
            %-------------------------------------------------------
            
            %         use only the untransformed image to get image power.
            mImg=softMaskPix.*img(pix);  % untransformed image
            I2=mImg(:)'*mImg(:);    % image power
            Ra2x=repmat(Ra2',ntai,1);  % ref power
            
            diffIR=I2-2*imgAmp*ccIR+imgAmp^2*Ra2x;  % squared difference between
            %         img and ref, ntAI=ntc*nAlphas*nIsos, nRefs*nVols
            
            if ri.usePrior
                % Create the WC log priors, which will be ntAI x nRefs in
                % size.
                logProbWC=reComputeLogProbWC5(nt,vesR,y0,si.mbnOffset,ri,moi,iImg);
            else
                logProbWC=0;
            end;
            wcPriors=reshape(logProbWC,ntai,narv)+repmat((logPVols(aVols)+logPRVs(arv(:)))',ntai,1);
                        
                if makeRawProbs
                    %         wcPriors has dimension(nTAI,nRefs,nVols)
                    wcPriors=reshape(wcPriors,ntc,nai,narv);
                    rawLogProbs.priors(tc,aip,arv(:),iImg)=wcPriors;
                    rawLogProbs.pXs(tc,aip,arv(:),iImg)=-reshape(diffIR,ntc,nai,narv)/(2*moi.sigmaN^2);
                end;
            wcPriors=reshape(wcPriors,ntai,narv);
            %   Compute the log of P(X|tarv), times the wcPriors
            logPXtairv=wcPriors-diffIR/(2*moi.sigmaN^2)-nPixEff*log(moi.sigmaN);
            %   exponentiate it and sum it to obtain logPX and logPtarv_X.
            maxLogPXtairv=max(logPXtairv(:));
            sclPXtairv=exp(logPXtairv-maxLogPXtairv);
            norm=sum(sclPXtairv(:));
            %         Get the full latent probabilities, p(t,a,i,r,v | X)
            ptairv_X=sclPXtairv/norm;
            logPX(iImg)=maxLogPXtairv+log(norm);  % log likelihood
            
            % get the separate probabilities
            ptai_X=sum(ptairv_X,2);  % sum over refs and vols
            %         pt_X=rotAlphaMat*ptai_X; % rotate and sum the translations
            pt_X0=sum(reshape(ptai_X,ntc,nai),2);  % no rotations
            pai_X=sum(reshape(ptai_X,ntc,nai),1)'; % p(alphaI): sum over trans
            pAlphas(aip,iImg)=pai_X; % p(alphaI): sum over trans
            prv_X=sum(ptairv_X,1)';  % sum over TAI to get p(r,v|X)
            pRefs(arv(:),iImg)=prv_X;
%             %         Get the best matches
            [~, ind]=max(ptairv_X(:));  % get 3D maximum, then get indices
%             for translations, alphaI for best match
            [indt, indai, indrv]=ind2sub([ntc nai nRV],ind);
            [tx,ty]=ind2sub([nt nt],tp(indt));
            bmTAI(:,iImg)=[tx ty aip(indai)]';
            bmClass(iImg)=rvp(indrv);
            % Get image expectation over trans and alpha. The weighting
            %         probability is ptairv_X which is ntai x nRV
            if useAltImgs
                eImg_tai_rv=transAltImgs*ptairv_X; % use altImgs for reconstruction
            else
                eImg_tai_rv=transImgs*ptairv_X; % result is nPix x nRV
            end;
            
            %             Check for out-of-bounds particle, in which case we don't
            %             accumulate its means
            if ~(any(tx==badTx) || any(ty==badTy) || any(aip(indai)==badAI(:)))
                % Accumulate the class means
                pixMeans(:,arv(:))=pixMeans(:,arv(:))+eImg_tai_rv*imgAmp;
                rawNorms(arv(:))=rawNorms(arv(:))+imgAmp^2*prv_X;  % nRV x 1
            else
                inBoundImgs(iImg)=false;
            end;
            
            %         Accumulate other model parameters
            localAmps=ccIR./Ra2x;
            imgAmps(iImg)=localAmps(:)'*ptairv_X(:);
            var=diffIR/nPixEff;
            eVar(iImg)=var(:)'*ptairv_X(:);  % expectation of variance
            pT_X0s(tp,iImg)=pt_X0;   % prob(translations)
            
        end; % if an active image
    end; % for iImg
    
    %     Create the full-size, ctf-filtered class means
    for k=1:nRV
        cMean=zeros(nw,nw);
        cMean(pix)=pixMeans(:,k);
        fclassMeans(:,:,k)=fclassMeans(:,:,k)+fft2(cMean).*ctfo;
    end;
    classNorms=classNorms+(ctf(:).^2)*rawNorms'; % nw^2 x nRV
end; % for iMi

% ------------ Output variables -----------
roi=struct;  % run output structure
roi.classMeans=reshape(real(ifft2(fclassMeans)),nw,nw,nRefs,nVols);
roi.classNorms=reshape(classNorms,nw,nw,nRefs,nVols);

%     alternate output variables
if isfield(ri.flags,'makeRawMeans') && ri.flags.makeRawMeans
    pixMeans=reshape(pixMeans,nPix,nRefs,nVols,nMi);
    rawMeans=zeros(nw^2,nRefs,nVols,nMi,'single');
    rawMeans(pix,:,:,:)=pixMeans;
    roi.rawMeans=reshape(rawMeans,nw,nw,nRefs,nVols,nMi);
    roi.rawNorms=reshape(rawNorms,nRefs,nVols,nMi);
end;
roi.intensiveFields=fieldnames(roi);  % All of the above are intensive

% Variables for each image (extensive variables)
roi.imgAmps=imgAmps;        % 1 x nImgs
roi.pTrans=pT_X0s;          % nt^2 x nImgs
roi.varNs=single(eVar);     % 1 x nImgs
roi.logPX=logPX;            % 1 x nImgs
roi.pRefs=reshape(pRefs,nGammas,nBetas,nVols,nImgs);
roi.pVols=shiftdim(sum(sum(roi.pRefs,1),2),2);            % nVols x nImgs
roi.pAlphas=reshape(pAlphas,nAlphas,nIsos,nImgs);       % nA x nI x nImgs
roi.inBoundImgs=inBoundImgs; % 1 x nImgs logical

% Best-match parameters
roi.imgClass=bmClass;  % best match ref for each image, nRV x nImgs
roi.imgTA=bmTAI; % (nt x nt x nAI x nim) best match trans and alpha

% alternate output variables
if isfield(ri.flags,'makeImgBestMatch') && ri.flags.makeImgBestMatch
    imgBestMatch=zeros(nw^2,nImgs,'single');
    imgBestMatch(pix,:)=pixBestMatch;
    roi.imgBestMatch=reshape(imgBestMatch,nw,nw,nImgs);
end;
if makeRawProbs
    roi.rawLogProbs=rawLogProbs;
end;
