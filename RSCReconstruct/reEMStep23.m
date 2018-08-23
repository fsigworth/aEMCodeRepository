function [classMeans,classNorms,roi,imgBestMatch,rawMeans,rawNorms,rawLogProbs]...
    =reEMStep23(imgs,refs,si,ri,moi,altImgs)
% function [classMeans, classNorms, roi,...]=reEMStep21(imgs,refs,si,ri,moi,altImgs)
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
% moi.activeTrans is nt x nt x nImgs booleans indicating which translations will be
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
useAltImgs=(nargin>5) && (numel(altImgs)>0) && ri.useAltImgs;
makeRawProbs=nargout>6;

sz=size(imgs);
nw=sz(1);     % Working image size
preRotated=numel(sz)>3;
nImgs=sz(end);

nt=ri.nTrans;    % Size of translations calculation, must be odd
nHalfTrans=floor((nt-1)/2);
nt2=nt^2;               % total trans. pixels to calculate
nT=nt2;

nMi=size(si.ctfs,3);    % number of micrographs

nRefs=size(refs,3);
nVols=max(1,numel(moi.pVols));

alphas=ri.alphas;  % we now give a list of alpha values, not including iso flips.
nAlphas=numel(ri.alphas);
nIsos=numel(ri.isos);
alphasI=ri.alphasI;

nRV=nRefs*nVols;
nAI=nAlphas*nIsos;
% % nTAI=nt2*nAI;

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
pAlphasI=zeros(nAI,nImgs,'single');
pRefs=zeros(nRV,nImgs,'single');

% Allocate the other output variables
logPX=zeros(1,nImgs);  % Likelihood.  Note that unassigned values are zero!
imgAmps=zeros(1,nImgs,'single');
eVar=zeros(1,nImgs,'single');  % noise variance
pVols=zeros(nVols,nImgs);
pT_X0s=zeros(nt2,nImgs,'single');  % prob of (absolute) translations
bmClass=zeros(1,nImgs,'single');  % best-matching ref/volume
bmTAI=zeros(3,nImgs,'single'); % best matching trans/alpha/iso
pixBestMatch=zeros(nPix,nImgs,'single'); % best transformed image

pixMeans=zeros(nPix,nRV,nMi);
rawNorms=zeros(nRV,nMi);
ctRefs=zeros(nPix,nRefs,nVols,nMi);
R2=zeros(nRefs,nVols,nMi);  % 2nd moment of refs
%********
if makeRawProbs
    rawLogProbs.priors=zeros(nT,nAI,nRV,nImgs,'single');
    rawLogProbs.pXs=zeros(nT,nAI,nRV,nImgs,'single');
end;
%*********

% Enumerate all the alphas including inside-out flips
% alphasI=zeros(nAlphas,nIsos,'single');
% for isi=1:nIsos
%     io=ri.isos(isi);
%     alphasI(:,isi)=alphas(:)+180*io;
% end;
% alphasI=alphasI(:)';

%     create the rotation matrix for translation
% rotAlphaMat=reMakeTransAlphaMatrix(alphasI,nt);

%         Create the ctf-filtered references
for iMi=1:nMi
    ctf=si.ctfs(:,:,iMi);
    ctfo=ifftshift(ctf);
    %    create the references, filtered by ctfs.
    %     They are multiplied by softMask and sampled into pixels
    for iVol=1:nVols;
        for iRef=1:nRefs;
            fref=fftn(refs(:,:,iRef,iVol));
            ctRef0=real(ifftn(fref.*ctfo));
            dc=invSoftMask(:)'*ctRef0(:);
            ctRef=ri.softMask.*(ctRef0-dc);
            ctRefPix=ctRef(pix);
            R2(iRef,iVol,iMi)=ctRefPix'*ctRefPix;
            ctRefs(:,iRef,iVol,iMi)=ctRefPix;
        end;
    end;
end;
ctRefs=reshape(ctRefs,nPix,nRV,nMi);
R2=reshape(R2,nRV,nMi);

% ---------Main loop over images---------
% % parfor iImg=1:nImgs  % actually runs slower!
for iImg=1:nImgs
    imgAmp=moi.imgAmps(iImg);
    tc=moi.activeTrans(:,:,iImg);
    if imgAmp>0 && any(tc(:))  % This image is active
        
        % Make the pointers to sparse translations
        tc=tc(:);  % active flags: booleans
        tp=find(tc);  % pointers of the active -> all translations
        ntc=numel(tp);  % number of active translations
        ntAI=ntc*nAI;
        
        vesR=si.rVesicle(iImg);
        y0=si.yClick(iImg);
        iMi=si.miIndex(iImg);  % micrograph (ctf) index

        %         ------set up the images------
        %   Images are rotated ccw.  That means if particle is rotated to
        %   right (cw rotated) then a positive alpha is the one that brings the particle into
        %   register.
        if preRotated
            unrotImgs=imgs(:,:,:,iImg);
            %           Get the power from an image with angle near 0
            avals=mod(alphas,90);
            avals(avals>45)=avals(avals>45)-90;
            [~,iAlpha0]=min(abs(avals));
            img=unrotImgs(:,:,iAlpha0); % pick this un-rotated one for power calc.
        else
            img=imgs(:,:,iImg);
            unrotImgs=rsRotateImage(img,-alphasI);
        end;        
        %         Create the matrix of shifted and rotated images
        transImgs=zeros(nPix,ntc,nAI);  % 1D images nTrans*nAlpha
        %         Make the translated images, multipled by soft mask
        for i=1:nAI  % loop over alphas
            uImg=unrotImgs(:,:,i);
            uImg=uImg(:);
            for j=1:ntc  % loop over translations
                transImgs(:,j,i)=uImg(transPix(:,tp(j))).*softMaskPix;
            end;
        end;
        transImgs=reshape(transImgs,nPix,ntAI);  % trans change most rapidly.

%         ------set up for altImgs------
        if useAltImgs  % in case we're using alternative images for recon
            if preRotated
                unrotAltImgs=altImgs(:,:,:,iImg);
            else
                unrotAltImgs=rsRotateImage(altImgs(:,:,iImg),-alphasI);
            end;
            transAltImgs=zeros(nPix,ntc,nAI);
            for i=1:nAI  % loop over alphas
                uAltImg=unrotAltImgs(:,:,i);
                uAltImg=uAltImg(:);
                for j=1:ntc  % loop over translations
                    transAltImgs(:,j,i)=uAltImg(transPix(:,tp(j))).*softMaskPix;
                end;
            end;
            transAltImgs=reshape(transAltImgs,nPix,ntAI);  % trans change most rapidly.
        end;
        
         % -------- Here's the big matrix multiplication --------
        ccIR=transImgs'*ctRefs(:,:,iMi);  % cross-correlation, ntAI x nRV
        %-------------------------------------------------------
        
        %         use only the untransformed image to get image power.
        mImg=softMaskPix.*img(pix);  % untransformed image
        I2=mImg(:)'*mImg(:);    % image power
        R2x=repmat(R2(:,iMi)',ntAI,1);  % ref power
        
        diffIR=I2-2*imgAmp*ccIR+imgAmp^2*R2x;  % squared difference between 
%         img and ref, ntAI=ntc*nAlphas*nIsos, nRefs*nVols
        
        % Create the WC log priors, which will be ntAI x nRefs in
        % size.
        logProbWC=reComputeLogProbWC3(nt,vesR,y0,si.mbnOffset,ri,moi,iImg);
        
        wcPriors=zeros(ntAI,nRefs,nVols);
        for iVol=1:nVols     % expand width by nVols
            if numel(moi.pRefs)>1
                wcPriors(:,:,iVol)=reshape(logProbWC,ntAI,nRefs)...
                +log(moi.pVols(iVol)+repmat(moi.pRefs(iImg,:,iVol),ntAI,1,1));
            else  % is this right??
                wcPriors(:,:,iVol)=reshape(logProbWC,ntAI,nRefs)...
                +log(moi.pVols(iVol)+moi.pRefs);
            end;
        end;
        if makeRawProbs
%         wcPriors has dimension(nTAI,nRefs,nVols)
            wcPriors=reshape(wcPriors,ntc,nAI,nRV);
            rawLogProbs.priors(tc,:,:,iImg)=wcPriors;
            rawLogProbs.pXs(tc,:,:,iImg)=-diffIR/(2*moi.sigmaN^2);
        end;
        wcPriors=reshape(wcPriors,ntAI,nRV);
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
        pt_X0=sum(reshape(ptai_X,ntc,nAI),2);  % no rotations
        pAlphasI(:,iImg)=sum(reshape(ptai_X,ntc,nAI),1)'; % p(alphaI): sum over trans
        prv_X=sum(ptairv_X,1);  % sum over TAI to get p(r,v|X)
        pRefs(:,iImg)=prv_X;
        %         [maxPrv, bmClass(iImg)]=max(prv_X);
        pVols(:,iImg)=sum(reshape(prv_X,nRefs,nVols),1);
%         Get the best matches
        [~, ind]=max(ptairv_X(:));  % get 3D maximum, then convert translations
        [indt, indAI, bmClass(iImg)]=ind2sub([ntc nAI nRV],ind);
        [tx,ty]=ind2sub([nt nt],tp(indt));
        bmTAI(:,iImg)=[tx ty indAI]';
        % Get image expectation over trans and alpha. The weighting 
%         probability is ptairv_X is nTAI x nRV
        if useAltImgs
            eImg_tai_rv=transAltImgs*ptairv_X; % use altImgs for reconstruction
        else
            eImg_tai_rv=transImgs*ptairv_X; % result is nPix x nRV
        end;
        % Accumulate the class means
        pixMeans(:,:,iMi)=pixMeans(:,:,iMi)+eImg_tai_rv*imgAmp;
        rawNorms(:,iMi)=rawNorms(:,iMi)+imgAmp^2*prv_X(:);  % nMi x nRV
%         Accumulate other model parameters
        localAmps=ccIR./R2x;
        imgAmps(iImg)=localAmps(:)'*ptairv_X(:);
        var=diffIR/nPixEff;
        eVar(iImg)=var(:)'*ptairv_X(:);  % expectation of variance
        pT_X0s(tp,iImg)=pt_X0;   % prob(translations)
        
    end; % if
end; % for iImg

%     Accumulate the classMeans and classNorms using the ctfs
fclassMeans=zeros(nw,nw,nRV,'single');  % will be complex
classNorms=zeros(nw*nw,nRV,'single');

for iMi=1:nMi
    ctf=si.ctfs(:,:,iMi);  % ctfs are zero-frequency in the center
    ctfo=ifftshift(ctf);
%     Create the full-size, ctf-filtered class means
    for k=1:nRV
        cMean=zeros(nw,nw);
        cMean(pix)=pixMeans(:,k,iMi);
        fclassMeans(:,:,k)=fclassMeans(:,:,k)+fft2(cMean).*ctfo;
    end;
    classNorms=classNorms+kron(ctf(:).^2,rawNorms(:,iMi)'); % nw^2 x nRV
end;
classMeans=reshape(real(ifft2(fclassMeans)),nw,nw,nRefs,nVols);
classNorms=reshape(classNorms,nw,nw,nRefs,nVols);

% Variables for each image
roi.imgAmps=imgAmps;        % 1 x nImgs
roi.pTrans=reshape(pT_X0s,nt,nt,nImgs);
roi.varNs=single(eVar);     % 1 x nImgs
roi.pVols=pVols;            % nVols x nImgs
roi.logPX=logPX;            % 1 x nImgs
roi.pRefs=reshape(pRefs,nRefs,nVols,nImgs);
roi.pAlphas=pAlphasI;       % nAI x nImgs

% Best-match parameters
roi.imgClass=bmClass;  % best match ref for each image, nRV x nImgs
roi.imgTA=bmTAI; % (nt x nt x nAI x nim) best match trans and alpha

%     other output variables
if nargout>4
    pixMeans=reshape(pixMeans,nPix,nRefs,nVols,nMi);
    rawMeans=zeros(nw^2,nRefs,nVols,nMi,'single');
    rawMeans(pix,:,:,:)=pixMeans;
    rawMeans=reshape(rawMeans,nw,nw,nRefs,nVols,nMi);
    rawNorms=reshape(rawNorms,nRefs,nVols,nMi);
end;
if nargout>3
    imgBestMatch=zeros(nw^2,nImgs);
    imgBestMatch(pix,:)=pixBestMatch;
    imgBestMatch=reshape(imgBestMatch,nw,nw,nImgs);
end;

