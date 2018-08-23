function [classMeans, classNorms, roi, pTrans, pAlphas, pRefs]=reEMStep21(imgs,refs,si,ri)
% function [classMeans, classNorms, roi, pTrans, pAlphas pRefs]=reEMStep21(imgs,refs,si,ri)

% new field: ri.imgActive, ri.radius (upper bound on hard mask radius), ri.softMask, ri.ctfs
% (all are assumed to have the same size, nw, as imgs and refs)

nw=size(imgs,1);     % Working image size
nt=ri.nTrans;    % Size of translations calculation
nHalfTrans=floor((nt-1)/2);
nt0=2*nHalfTrans+1;     % actual no. of points to calculate, an odd number
nt2=nt0^2;               % total trans. pixels to calculate

nImgs=size(imgs,3);
nRefs=size(refs,3);
nVols=numel(ri.pVols);
nRV=nRefs*nVols;

nAlphas=ri.angleN(1);
nBetas=ri.angleN(2);
nGammas=ri.angleN(3);
alphas=ri.angleMin(1)+ri.angleStep(1)*(0:nAlphas-1);


% Create the hard mask for selecting pixels
% pick the radius for the hard mask to be maximum to allow translations
hmRadius=floor(min(ri.radius,floor(nw-1)/2-nHalfTrans));


% Mask for centered particles
msk=fuzzymask(nw,2,hmRadius,0); % should be integer, 1 more than radius
invSoftMask=1-ri.softMask;
invSoftMask=invSoftMask/sum(invSoftMask(:)); % sums to 1
nPixEff=ri.softMask(:)'*ri.softMask(:); % softMask has pixels = 1

pix=msk(:);  % nw^2 x 1: picks out the pixels for centered mask
nPix=sum(pix);
softMaskPix=ri.softMask(pix);

% For particles with all shifts
transPix=false(nw^2,nt0^2);
offset=floor(nw/2+1)-nHalfTrans-1;
for ix=1:nt0
    for iy=1:nt0
        msk=fuzzymask(nw,2,hmRadius,0,[ix iy]+offset);
        transPix(:,ix+(iy-1)*nt0)=msk(:);
    end;
end;
% Check for consistency
if ~all(sum(transPix,1)==sum(pix))
    error('Hard mask size varies with translations');
end;

% Allocate the output variables
pAlphas=zeros(nAlphas,nImgs,'single');
pRefs=zeros(nRV,nImgs,'single');
% pta_Xs=zeros(nt2,nAlphas,nImgs,'single');

% Allocate the other output variables
logPX=-inf*ones(nImgs,1);
% pTrans=zeros(nt0,nt0,nImgs,'single');
imgAmps=zeros(nImgs,1);
% eAmps=zeros(nImgs,1);
eVar=zeros(nImgs,1);
pVols=zeros(nVols,nImgs);
pt_Xs=zeros(nt2,nImgs);

fclassMeans=zeros(nw,nw,nRV,'single');  % will be complex
classNorms=zeros(nw*nw,nRV,'single');



% Loop over micrographs (=defocus groups)
nMi=size(si.ctfs,3);
for iMi=1:nMi
    iImgs=find(si.miIndex==iMi);
    %     nImgs_ctf=numel(iImgs);
    ctf=si.ctfs(:,:,iMi);
    ctfo=ifftshift(ctf);
    %     accumulators for 1 micrograph
    eImg_ta_rvu=zeros(nPix,nRefs*nVols);  % accumulate class means in defocus groups
    eNorm_rvu=zeros(nRefs*nVols,1);
    
    %    create the references
    ctRefs=zeros(nPix,nRefs,nVols);
    R2=zeros(nRefs,nVols);  % 2nd moment of refs
    for iVol=1:nVols;
        for iRef=1:nRefs;
            fref=fftn(refs(:,:,iRef,iVol));
            ctRef0=real(ifftn(fref.*ctfo));
            dc=invSoftMask(:)'*ctRef0(:);
            ctRef=ri.softMask.*(ctRef0-dc);
            ctRefPix=ctRef(pix);
            R2(iRef,iVol)=ctRefPix(:)'*ctRefPix;
            ctRefs(:,iRef,iVol)=ctRefPix;
        end;
    end;
    ctRefs=reshape(ctRefs,nPix,nRV);
    
    %     create the rotation matrix for translation
    rotAlphaMat=reMakeTransAlphaMatrix(alphas,nt0);
    
    for iImg=iImgs'  % loop over all particles from this micrograph
        
        imgAmp=ri.imgAmps(iImg);
        if imgAmp>0 && ri.imgActive(iImg) % This image is active
            img0=imgs(:,:,iImg);
            imgDC=invSoftMask(:)'*img0(:);
            img=ri.softMask.*(img0-imgDC);
            vesR=si.rVesicle(iImg);
            y0=si.yClick(iImg);
            
            unrotImgs=rsRotateImage(img,-alphas);
            %         Create the matrix of shifted and rotated images
            transImgs=zeros(nPix,nt2,nAlphas);  % 1D images nTrans*nAlpha
            for i=1:nAlphas
                uImg=ri.softMask.*unrotImgs(:,:,i);  % softmask is not shifted, ok??*****
                uImg=uImg(:);
                for j=1:nt2
                    transImg=uImg(transPix(:,j)).*softMaskPix;
                    %                     imgPower(iImg)=transImg'*transImg;
                    transImgs(:,j,i)=transImg;
                end;
            end;
            transImgs=reshape(transImgs,nPix,nt2*nAlphas);  % trans change most rapidly.
            
            % -------- Here's the big matrix multiplication --------
            ccIR=transImgs'*ctRefs*imgAmp;  % cross-correlation, nt2*nAlphas x nRefs*nVols
            
            imgPower=sum(transImgs.^2,1)';
            R2x=repmat(imgAmp^2*R2(:)',nt2*nAlphas,1);
            I2x=repmat(imgPower,1,nRefs*nVols);
            
            diffIR=I2x-2*ccIR+R2x;  % difference between img and ref, nt2*nAlphas, nRefs*nVols
            
            % Create the WC log priors, which will be nt2*nAlphas x nRefs in
            % size.
            wcPriors1vol=zeros(nt2*nAlphas,nRefs);
            betaStep=ri.angleStep(2);
            for i=1:nBetas
                beta=ri.angleMin(2)+(i-1)*betaStep;
                dBeta=4*pi*betaStep/180*sind(beta);
                logPwc=reComputeLogProbWC(nt0,ri.sigmaC,ri.sigmaG,vesR,...
                    y0,alphas,1/nAlphas,beta,dBeta);
                %             nt0 x nt0 x nAlphas in size
                betaOffset=(i-1)*nGammas;
                wcPriors1vol(:,1+betaOffset:nGammas+betaOffset)...
                    =repmat(logPwc(:),1,nGammas);
            end;
            wcPriors=zeros(nt2*nAlphas,nRefs,nVols);
            for i=1:nVols     % expand width by nVols
                wcPriors(:,:,iVol)=wcPriors1vol+log(ri.pVols(iVol));
            end;
            wcPriors=reshape(wcPriors,nt2*nAlphas,nRefs*nVols);
            
            %             Compute the log of P(X|tarv), times the wcPriors
            logPXtarv=wcPriors-diffIR/(2*ri.sigmaN^2)-nPixEff*log(ri.sigmaN);
            % exponentiate it and sum it to obtain logPX and logPtarv_X.
            maxLogPXtarv=max(logPXtarv(:));
            sclPXtarv=exp(logPXtarv-maxLogPXtarv);
            norm=sum(sclPXtarv(:));
            ptarv_X=sclPXtarv/norm;        % full latent probabilities
            logPX(iImg)=maxLogPXtarv+log(norm);  % log likelihood
            
            % get the separate probabilities
            pta_X=sum(ptarv_X,2);  % sum over refs and vols
            pt_X=rotAlphaMat*pta_X; % rotate and sum the translations
            pAlphas(:,iImg)=sum(reshape(pta_X,nt2,nAlphas),1)'; % p(alpha): sum over trans
            prv_X=sum(ptarv_X,1);  % sum over trans and alpha
            pRefs(:,iImg)=prv_X;
            pVols(:,iImg)=sum(reshape(prv_X,nRefs,nVols),1);
            % get expectation values
            %         transImgs is nPix x nt2*nImgAlphas
            %         ptarv_X is nt2*nImgAlpha x nRefs x nVols
            %             Accumulate the class means
            eImg_ta_rv=transImgs*ptarv_X; % result is nPix x nRefs*nVols
            cMeans=zeros(nw^2,nRV,'single');
            cMeans(pix,:)=eImg_ta_rv;
            ctfos=repmat(ctfo,1,1,nRV);
            fclassMeans=fclassMeans+fft2(reshape(cMeans,nw,nw,nRV)).*ctfos*imgAmp;
            classNorms=classNorms+kron(ctf(:).^2*imgAmp^2,prv_X(:)'); % nw^2 x nRV
            
%             eAmps(iImg)=imgAmps(:)'*ptarv_X(:);
            localAmps=ccIR./R2x;
            imgAmps(iImg)=localAmps(:)'*ptarv_X(:); % need to divide by imgAmp? ********
            var=diffIR.^2;
            eVar(iImg)=var(:)'*ptarv_X(:);
            
            %             Accumulate the class means
            %             eImg_ta_rvu=eImg_ta_rvu+eImg_ta_rv;
            %             eNorm_rvu=eNorm_rvu+prv_X(:);  % nRefs*nVols x 1
            
            %             pta_Xs(:,iImg)=pta_X;
            pt_Xs(:,iImg)=pt_X;
            
        end; % if
    end; % for iImg
    %     Operate with CTF on classes for the defocus group
    %     cMeans=zeros(nw,nw,nRef*nVols,'single');
    %     cMeans(pix,:)=eImg_ta_rvu;
    %     cMeans=reshape(cMeans,nw,nw,nRefs*nVols);
    %     classMeans=classMeans+real(ifft2(fft2(cMeans).*ctfos));
    %     classNorms=classNorms+kron(ctf(:),eNorm_rvu'); %  nw^2 x nRefs*nVols
end; % parfor iMi
classMeans=reshape(real(ifft2(fclassMeans)),nw,nw,nRefs,nVols);
classNorms=reshape(classNorms,nw,nw,nRefs,nVols);
roi.imgAmps=imgAmps;
pTrans=reshape(pt_Xs,nt,nt,nImgs);
roi.sigmaNs=eVar;
roi.pVols=pVols;
pRefs=reshape(pRefs,nRefs,nVols,nImgs);