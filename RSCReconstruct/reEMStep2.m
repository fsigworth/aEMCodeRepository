function [classMeans, classNorms, roi, pTrans, pAlphas pRefs]=reEMStep2(imgs,refs,emFlags,si,ri,moi)

% function [cMeans,cNorms,roi,pTrans,pAngles]=reEMStep2(imgs,refs,emFlags,ri,moi)
%  Abstraction of the main EM step for single-particle reconstruction.
% --Inputs are stacks of images and references, along with the very large
% boolean array emFlags (nAlphas x nRefs x nVols x nImgs) which tells
% which alpha angles and references to use for each image. We use the
% run-info structure ri fields
%   ri.angleMin, ri.angleStep, ri.angleN (for computing angles)
%   and ri.n sets the working image size, ri.nt is the size of translations
%   to compute.
% We use the si fields si.ctfs, si.miIndex, si.rVesicle, si.yClick.
% We use the model info structure fields moi. sigmaG, moi.sigmaC,
%   moi.sigmaN, moi.pVols, moi.imgAmps
%
% --Outputs are the classMeans (real images) and classNorms (zero-center
% Fourier) Along with the roi (run-info) outputs for each image,
%   roi.logLik, roi.pVols, roi.imgAmps
% and the parameter expectation values ePars.  pTrans is an nt x nt x nImgs
% stack of translation probabilities, pAlphas is the probability of alpha
% angles and pRefs is a matrix giving ref indices for each image. Its size
% is nRefs*nVols x nImgs

% ri.n, ri.radius, ri.nTrans


% Variable naming examples
% pta_rvX = Prob(trans, alpha | ref, vol, X),
%   has dimension nw x nw x nAlpha x nRef x nVol.
% ePars_a_rv = Expectation_a{Pars}(ref,vol),
%   has dimension nPars x nRef x nVol.
% There are subsets of the alpha indices for alpha values that we actually
% use.  imgAlphaInds are ones that are used for a given image.
% refAlphaFlags indicates which of these are actually used for a given
% image and ref combination.  In naming variables we use a1 and a2 to
% indicate these two kinds of pointers, for example the line
%  pa1_rvX(refAlphaFlags,k,iVol)=pa2_rX; % nImgAlphas x nRefs in size


nPars=5;

nw=ri.n;                % Working image size
nt=ri.nTrans;    % Size of translations calculation
nt2=nt^2;
nHalfTrans=floor((nt-1)/2);
nt0=2*nHalfTrans+1;     % actual no. of points to calculate

ds=size(imgs,1)/nw;     % Downsampling factor

[nAlphas,nImgs]=size(alphaFlags);
[nRefs,nVols,nImgs]=size(refFlags);


% Create the hard mask for selecting pixels
% pick the radius for the hard mask to be maximum
radius=floor(min(ri.radius,floor(nw-1)/2-nHalfTrans));
% Mask for centered particles
msk=fuzzymask(nw,radius,0);
pix=msk(:);  % nw^2 x 1: picks out the pixels for centered mask
nPix=sum(pix);
% For particles with all shifts
translatedPix=false(nw^2,nt0^2);
offset=floor((nw-1)/2)-nHalfTrans;
for ix=1:nt0
    for iy=1:nt0
        msk=fuzzymask(nw,radius,0,[ix iy]+offset);
        transPix(:,ix+(iy-1)*nt0)=msk(:);
    end;
end;
% Check for consistency
if ~all(sum(transPix,1)==sum(pix))
    error('Hard mask size varies with translations');
end;

% Create the references
ctRefs=zeros(nPix,nRefs,nVols,nCTFs);
for iVol=1:nVols;
    
    for iRef=1:nRefs;
        fref=Cropo(fftn(refs(:,:,iRef,iVol)),nw);
        for i=1:nCTFs
            ctRef=ri.softMask.*real(ifftn(fref.*ifftshift(Crop(ctf,nw))));
            ctRefs(:,iRef,iVol,i)=ctRef(pix);
        end;
    end;
end;

% Allocate the big output variables
pAlphas=zeros(nAlphas,nImgs,'single');
pRefs=zeros(nRefs,nVols,nImgs,'single');  % someday make this uint8?

% Allocate the other output variables
logLik=-inf*ones(nImgs,1);
ePars=zeros(nImgs,nPars);
pTrans=zeros(nt,nt,nImgs,'single');
imgAmps=zeros(nImgs,1);
pVols=zeros(nImgs,nVols);
classMeans=zeros(nw,nw,nRefs,nVols);
classNorms=zeros(nw,nw,nRefs,nVols);

eImg_rvs=zeros(nPix,nRefs*nVols,nCTFs);  % accumulate class means in defocus groups
eNorm_rvs=zeros(nRefs,nVols,nCTFs);


% for iImg
%     
%    downsample the image, rotate it, maybe pad it. 
%    if new ctf, filter the refs and convert to 1D
%    convert rotated and transl image to 1D
%    form transformedImage x refs = nTransformed x nRefs*nVols matrix
%    copy and divide by R2s to get est. ampl.
%    add I2 and R2 to get DiffIR
%    add priors to this matrix to get logp(tarv_X)
%    Exponentiate to convert to p(tarv_X) and likelihood p(X).
%    compute p(tarv_X)(:)*[DiffIR(:)' EstAmpl(:)'] to get expect. vals.
%    Sum p(tarv_X) over nTransformed to get p(rv_X)
%    Sum p(tarv_X) over nRefs*nVols to get p(ta_X)
%    Sum p(ta) over trans to get p(alpha_X)
%    Rotate and sum bands of p(ta) to get p(t_X) relative to original image.
%    

parfor iImg=1:nImgs
    %     Downsample the image and variables to the working size nw
    imgAmp=ri.imgAmps(iImg);
    imgRefFlags=refFlags(:,:,iImg);% logical(nalphas x nrefs x nvols)
    imgRefInds=find(imgRefFlags);
    imgAlphaFlags=alphaFlags(:,iImg);
    imgAlphaInds=find(imgAlphaFlags);  % all alphas used by this image
    nImgAlphas=numel(imgAlphaInds);
    if imgAmp>0 && numel(imgAlphaInds)>0 && any(imgRefFlags(:)) % This image is active
        nPixEff=ri.softMask(:)*ri.softMask(:)';
        img=Downsample(imgs(:,:,iImg),nw);
        ctfIndex=si.miIndex(iImg);
        vesR=si.rVesicle(iImg)/ds;
        y0=si.yClick(iImg)/ds;
        imgAlphas=ri.angleMin(1)+(imgAlphaInds-1)*ri.angleStep(1); % actual angles
        unrotImgs=Crop(rsRotateImage(img,-imgAlphas)),nw,1;
        transImgs=zeros(nPix,nt2,nImgAlphas);
        for i=1:nImgAphas
            uImg=ri.softMask.*unrotImgs(:,:,i);
            uImg=uImg(:);
            for j=1:nt2
                transImgs(:,j,i)=uImg(transPix(:,j));
            end;
        end;
        transImgs=reshape(transImgs,nPix,nt2*nImgAlphas);  % trans change most rapidly.
        refArray=ctRefs(:,imgRefFlags,ctfIndex);

      % Here's the big matrix multiplication
        ccIR=transImgs*refArray;  % cross-correlation nImgAlphas*n2 x nRefs*nVols
        
        R2=zeros(nImgRefs,1);
        for i=1:nImgRefs
            R2(i)=refArray(:,i)*refArray(:,i)';
        end;
        I2=zeros(nImgAlphas*nt2,1);
        for i=1:nImgAlphas*nt2
            I2(i)=transImgs(:,i)*transImgs(:,i)';
        end;
        R2x=repmat(R2',nImgAlphas*nt2,1);
        I2x=repmat(R1,1,nImgRefs*nVols);
        
        amps=ccIR./R2x;
        diffIR=I2x-2*ccIR+R2x;
        
%         Create the WC log priors
%         Get the distinct beta values
        wcPriors1vol=zeros(nt2*nImgAlphas,nImgRefs);
        betaInds=idivide(imgRefInds,ri.angleN(3));
        betaSteps=find([1 diff(betaInds) 1]~=0);
        for i=betaSteps
            beta=(betaInds(i)-1)*ri.angleStep(2);
            dBeta=4*pi*ri.angleStep(2)/180*sind(beta);
            logPwc=reComputeLogProbWC(nt0,ri.sigmaC,ri.sigmaG,si.vesicleR(iImg),...
                si.y0(iImg),imgAlphas,1/ri.angleN(1),beta,dBeta);
            wcPriors1vol(:,betaInds(i):betaInds(i+1)-1)...
                =repmat(logPwc(:),1,betaInds(i+1)-betaInds(i));
        end;
        wcPriors=zeros(nt2*nImgAlphas,nImgRefs,nVols);
        for i=1:nVols
        
            wcPriors(:,:,iVol)=wcPriors1vol+log(ri.pVol(iVol));
        end;
        wcPriors=reshape(wcPriors,nt2*nImgAlphas,nImgRefs*nVols);
        logPXtarv=wcPriors-diffIR/(2*ri.sigmaN^2)-nPixEff*log(sigmaN);
% exponentiate it and sum it to obtain logPX and logPtarv_X.
        maxLogPXtarv=max(logPXtarv(:));
        sclPXtarv=exp(logPXtarv-maxLogPXtarv);
        norm=sum(sclPXtarv(:));
        ptarv_X=sclPXtarv/norm;
        logPX=maxLogPXtarv+log(norm);
        
% get the separate probabilities
        pta_X=sum(ptarv_X,2);  % sum over refs and vols
        pa_X=sum(reshape(pta_X,nt2,nImgAlphas),1);
        prv_X=sum(ptarv_X,1);  % sum over trans and alpha
        
% get expectation values
%         transImgs is nPIx x nt2*nImgAlphas
%         ptarv_X is nt2*nImgAlpha x nRefs x nVols
        eImg_ta_rv=transImgs*ptarv_X; % result is nPix x nRefs*nVols
        eAmps(iImg)=amps(:)*ptarv_X(:);
        var=diffIR.^2;
        eVar(iImg)=var(:)*ptarv_X(:);

        eImg_ta_rvs(:,:,ctfIndex)=eImg_ta_rvs(:,:,ctfIndex)+eImg_ta_rv;
        eNorm_rvs(:,ctfIndex)=eNorm_rvs(:,ctfIndex)+prv_X;
        
        pta_Xs(:,iImg)=pta_X;  % need to index alphas
        
%         need to index refs
    end; % if
end; % parfor iImg

% unrotate the ptas
% sum and filter the class means

    
        
