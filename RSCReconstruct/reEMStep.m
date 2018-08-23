function [classMeans, classNorms, roi, pTrans, pAngles]=reEMStep(imgs,refs,emFlags,si,ri,moi)
% function [cMeans,cNorms,roi,pTrans,pAngles]=reEMStep(imgs,refs,emFlags,ri,moi)
%  Abstraction of the main EM step for single-particle reconstruction.
% --Inputs are stacks of images and references, along with the very large
% boolean structure emFlags (nAlphas x nRefs x nVols x nImgs) which tells
% which alpha angles and references to use for each image. We use the
% run-info structure ri fields
%   ri.angleMin, ri.angleStep, ri.angleN (for computing angles)
%   and ri.n sets the working image size
% We use the si fields si.ctfs, si.miIndex, si.rVesicle, si.yClick.
% We use the model info structure fields moi. sigmaG, moi.sigmaC,
%   moi.sigmaN, moi.pVols, moi.imgAmps
%
% --Outputs are the classMeans (real images) and classNorms (zero-center
% Fourier) Along with the roi (run-info) outputs for each image,
%   roi.logLik, roi.pVols, roi.imgAmps
% and the parameter expectation values ePars.  pTrans is an nw x nw x nImgs
% stack of translation probabilities, and pAngles is a sparse matrix giving
% all angle probabilities for each image. The matrix size is the same
% overall size and has the same zeros as emFlags. It has dimensions
% nalphas*nRefs*nVols x nImgs;

% Uses the functions
%   rsRotateImage
%   reComputeLogProbWC  % Geometry prior
%   reOneImgOneRef      % Inner loop

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
ds=size(imgs,1)/nw;     % Downsampling factor

[nAlphas,nRefs,nVols,nImgs]=size(emFlags);
% [nAlphas,nRefs,nVols,nImgs,is,js]=reGetFlagSizes(emFlags);

% Allocate the big output variable
pAngles=zeros(nAlphas,nRefs,nVols,nImgs,'single');  % someday make this uint8?

% Allocate the other output variables
logLik=-inf*ones(nImgs,1);
ePars=zeros(nImgs,nPars);
pTrans=zeros(nw,nw,nImgs,'single');
imgAmps=zeros(nImgs,1);
pVols=zeros(nImgs,nVols);
classMeans=zeros(nw,nw,nRefs,nVols);
classNorms=zeros(nw,nw,nRefs,nVols);

parfor iImg=1:nImgs
    %     Downsample the image and variables to the working size nw
    imgAmp=ri.imgAmps(iImg);
    imgFlags=emFlags(:,:,:,iImg);% logical(nalphas x nrefs x nvols)
    imgAlphaInds=find(any(any(imgFlags,2),3));  % all alphas used by this image
    if imgAmp>0 && numel(imgAlphaInds)>0  % This image is active
        img=Downsample(imgs(:,:,iImg),nw);
        ctf=Crop(si.ctfs(:,:,si.miIndex(iImg)),nw);
        ctfo=ifftshift(ctf);  % zero frequency is at origin
        vesR=si.rVesicle(iImg)/ds;
        y0=si.yClick(iImg)/ds;
        
        % Allocate the image accumulators
        eImg_a_rv=zeros(nw,nw,nRefs,nVols);
        ePars_a_rv=zeros(nPars,nRefs,nVols);
        logPXrv=-inf*ones(nRefs,nVols);
        pa1_rvX=zeros(nAlphas,nRefs,nVols);
        
        %     Pick up the flags telling us which alphas and refs to search
        %     imgFlags=reGetFlagsForImage(emFlags,iImg);% logical(nalphas x nrefs x nvols)
        imgAlphas=ri.angleMin(1)+(imgAlphaInds-1)*ri.angleStep(1); % actual angles
        nImgAlphas=numel(imgAlphas);  % indexed as a1
        fUnrotImgsA1=zeros(nw,nw,nImgAlphas,'single');
        fUnrotImgsA1(:,:,imgAlphaInds)=Crop(fft2(rsRotateImage(img,-imgAlphas)),nw,1);
        
        pt_a1rvX=zeros(nw,nw,nImgAlphas,nRefs,nVols);
        parv_X=zeros(nAlphas,nRefs,nVols);
        %         Pick the volumes that are in use
        activeVols=find(squeeze(any(any(imgFlags,2),1)))';
        logProbWC_a1b=[];
        for iVol=activeVols  % Loop over volumes
            refInds=find(any(imgFlags(:,:,iVol),1));  % flags in the set of refs for this volume.
            oldBeta=inf;
            
            for k=refInds  % Loop over references
                %             Get the subset of alphas, indexed a2, that we'll use with this reference
                refAlphaFlags=imgFlags(imgAlphaInds,k,iVol);
                beta=reGetRefAngles(k,ri);
                pRef=ri.angleStep(2)/360*sind(beta);  % prior for this reference
                %             Get the geometry priors
                if beta~=oldBeta
                    logProbWC_a1b=reComputeLogProbWC(nw,moi.sigmaC,moi.sigmaG,...
                        vesR,y0,imgAlphas,ri.angleStep(1),beta,pRef);
                    oldBeta=beta;
                end;
                
                %         ------    Do the single-reference calculation   -----
                %             We use only the alpha values pointed to by refAlphaFlags (a2)
                cref=real(ifftn( Cropo(fftn(refs(:,:,k,iVol)),nw).*ctfo*imgAmp ));
                %             cref is scaled by imgAmp here.
                [logPX, accumImg, ePars_a, pa2_rX, pt_a2rX]...
                    =reOneImgOneRef(fUnrotImgsA1(:,:,refAlphaFlags),...
                      logProbWC_a1b(:,:,refAlphaFlags),moi.sigmaN,cref);
                
                %              Store the results
                
                eImg_a_rv(:,:,k,iVol)=real(ifftn( fftn(accumImg).*ctfo*imgAmp ));  % nw x nw x nRefs x nVols
                logPXrv(k,iVol)=logPX+log(pRef*moi.pVols(iVol)); % =p(X|rv)*p(X)*p(v)
                ePars_a_rv(:,k,iVol)=ePars_a;
                pa1_rvX(refAlphaFlags,k,iVol)=pa2_rX; % p(alpha|X,ref)      nImgAlphas x nRefs
                pt_a1rvX(:,:,refAlphaFlags,k,iVol)=pt_a2rX;  % nw x nw x nalpha x nRefs
            end; % k loop over refs
        end; % ivol loop over vols
        %     Compute p(rv_X)
        maxLogPXrv=max(logPXrv(:));
        sclPXrv=exp(logPXrv-maxLogPXrv);
        normPXrv=sum(sclPXrv(:));  % sum over both r and v
        logLik(iImg)=log(normPXrv)+maxLogPXrv;  % Lik=sum_r,v{p(X,r,v)}
        prv_X=sclPXrv/normPXrv;  %  nRefs x nVols
        pta_X=reshape(pt_a1rvX,nw^2*nImgAlphas,nRefs*nVols)*prv_X(:); % nw^2 x nImgAlphas
        
        %     output variables
        %     Rotate the translations by alpha and combine
        pTrans(:,:,iImg)=sum(rsRotateImage(reshape(pta_X,nw,nw,nImgAlphas),imgAlphas),3);
        tPars=reshape(ePars_a_rv,nPars,nRefs*nVols)*prv_X(:);
        tPars(2)=tPars(2)*imgAmp;  % correct for scaling of ctf
        imgAmps(iImg)=tPars(2);
        ePars(iImg,:)=tPars;
        
        ctNorm=(imgAmp*ctf).^2;
        
        parv_X(imgAlphaInds,:,:)...
            =repmat(reshape(prv_X,1,nRefs,nVols),nImgAlphas,1,1).*pa1_rvX;
        %     pAngles(:,iImg)=parv_X(:);
        pAngles(:,:,:,iImg)=parv_X;
        imgPrv_X=repmat(shiftdim(prv_X,-2),nw,nw,1,1);  % nw^2 x nRefs x nVols
        classMeans=classMeans+eImg_a_rv.*imgPrv_X;
        classNorms=classNorms+repmat(ctNorm,1,1,nRefs,nVols).*imgPrv_X;
        pVols(iImg,:)=sum(prv_X,1);
    end; % if imgAmp
end; % for iImg
% classMeans=reshape(classMeans,nw,nw,nRefs,nVols);
% classNorms=reshape(classNorms,nw,nw,nRefs,nVols);
roi.logLik=logLik;
roi.pVols=pVols;
roi.ePars=ePars;
roi.imgAmps=imgAmps;

