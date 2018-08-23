function [classMeans, classNorms, roi, pTrans, pAlphas, pRefs]=reEMStep22p(imgs,refs,si,ri)
% function [classMeans, classNorms, roi, pTrans, pAlphas pRefs]=reEMStep22(imgs,refs,si,ri)
% p version uses all pixels in the images.
%
% new field: ri.imgActive, ri.radius (upper bound on hard mask radius), ri.softMask, ri.ctfs
% (all are assumed to have the same size, nw, as imgs and refs)

nw=size(imgs,1);     % Working image size
nt=ri.nTrans;    % Size of translations calculation
nHalfTrans=floor((nt-1)/2);
nt0=2*nHalfTrans+1;     % actual no. of points to calculate, an odd number
nt2=nt0^2;               % total trans. pixels to calculate

nMi=size(si.ctfs,3);

nImgs=size(imgs,3);
nRefs=size(refs,3);
nVols=numel(ri.pVols);

nAlphas=ri.angleN(1);
nBetas=ri.angleN(2);
nGammas=ri.angleN(3);
alphas=ri.angleMin(1)+ri.angleStep(1)*(0:nAlphas-1);

nRV=nRefs*nVols;
nTA=nt2*nAlphas;

% Create the hard mask for selecting pixels
% pick the radius for the hard mask to be maximum to allow translations
% % hmRadius=floor(min(ri.radius,floor(nw-1)/2-nHalfTrans));


% Mask for centered particles
% % msk=fuzzymask(nw,2,hmRadius,0); % r should be integer, 1 more than actual radius
% % invSoftMask=1-ri.softMask;
% % invSoftMask=invSoftMask/sum(invSoftMask(:)); % sums to 1
nPixEff=ri.softMask(:)'*ri.softMask(:); % softMask has pixels approaching 1

% % pix=msk(:);  % nw^2 x 1: picks out the pixels for centered mask
% % nPix=sum(pix);
nPix=nw^2;

% % softMaskPix=ri.softMask(:);

% % % For particles with all shifts
% % transPix=false(nw^2,nt0^2);
% % offset=floor(nw/2+1)-nHalfTrans-1;
% % for ix=1:nt0
% %     for iy=1:nt0
% %         msk=fuzzymask(nw,2,hmRadius,0,[ix iy]+offset);
% %         transPix(:,ix+(iy-1)*nt0)=msk(:);
% %     end;
% % end;
% % % Check for consistency
% % if ~all(sum(transPix,1)==sum(pix))
% %     error('Hard mask size varies with translations');
% % end;

% Allocate the output variables
pAlphas=zeros(nAlphas,nImgs,'single');
pRefs=zeros(nRV,nImgs,'single');
% pta_Xs=zeros(nt2,nAlphas,nImgs,'single');

% Allocate the other output variables
logPX=zeros(nImgs,1);  % note that unassigned values are zero!
% pTrans=zeros(nt0,nt0,nImgs,'single');
imgAmps=zeros(nImgs,1);
% eAmps=zeros(nImgs,1);
eVar=zeros(nImgs,1);
pVols=zeros(nVols,nImgs);
pt_Xs=zeros(nt2,nImgs);

cMeans=zeros(nPix,nRV,nMi);
cNorms=zeros(nRV,nMi);
ctRefs=zeros(nPix,nRefs,nVols,nMi);
R2=zeros(nRefs,nVols,nMi);  % 2nd moment of refs

%     create the rotation matrix for translation
rotAlphaMat=reMakeTransAlphaMatrix(alphas,nt0);

%         Create the ctf-filtered references
for iMi=1:nMi
    %     nImgs_ctf=numel(iImgs);
    ctf=si.ctfs(:,:,iMi);
    ctfo=ifftshift(ctf);
    
    %    create the references, filtered by ctfs.
    %     They are multiplied by softMask and sampled into pixels
    for iVol=1:nVols;
        for iRef=1:nRefs;
            fref=fftn(refs(:,:,iRef,iVol));
            ctRef0=real(ifftn(fref.*ctfo));
            %             dc=invSoftMask(:)'*ctRef0(:);
            dc=0;
            ctRef=ri.softMask.*(ctRef0-dc);
            ctRefPix=ctRef(:);
            R2(iRef,iVol,iMi)=ctRefPix'*ctRefPix;
            ctRefs(:,iRef,iVol,iMi)=ctRefPix;
        end;
    end;
end;
ctRefs=reshape(ctRefs,nPix,nRV,nMi);
R2=reshape(R2,nRV,nMi);
%     f1=@AccumCMeans;
%     f2=@AccumCNorms;
%

% ---------Main loop over images---------
% parfor iImg=1:nImgs
for iImg=1:nImgs
    
    localCMeans=zeros(nPix,nRV,nMi);
    localCNorms=zeros(nRV,nMi);
    imgAmp=ri.imgAmps(iImg);
    if imgAmp>0 && ri.imgActive(iImg) % This image is active
        img=imgs(:,:,iImg);
        %         img0=imgs(:,:,iImg);
        %         imgDC=invSoftMask(:)'*img0(:);
        %         img=ri.softMask.*(img0-imgDC);
        vesR=si.rVesicle(iImg);
        y0=si.yClick(iImg);
        iMi=si.miIndex(iImg);  % micrograph (ctf) index
        
        unrotImgs=rsRotateImage(img,-alphas);
        %         Create the matrix of shifted and rotated images
        % %         transImgs=zeros(nPix,nt2,nAlphas);  % 1D images nTrans*nAlpha
        transImgs=zeros(nPix,nt0,nt0,nAlphas,'single');
        %         Make the translated images, multipled by soft mask
        for iA=1:nAlphas
            % %             uImg=unrotImgs(:,:,i);
            % %             uImg=uImg(:);
            % %             for j=1:nt2
            % %                 transImg=uImg(transPix(:,j)).*softMaskPix;
            % %                 %                     imgPower(iImg)=transImg'*transImg;
            % %                 transImgs(:,j,i)=transImg;
            % %             end;
            um=unrotImgs(:,:,iA);
            for jy=1:nt0
                for jx=1:nt0
                    uImg=ri.softMask.*circshift(um,[jx jy]-nHalfTrans-1);
                    transImgs(:,jx,jy,iA)=uImg(:);
                end;
            end;
        end;
        transImgs=reshape(transImgs,nPix,nTA);  % trans change most rapidly.
        
        % -------- Here's the big matrix multiplication --------
        ccIR=transImgs'*ctRefs(:,:,iMi)*imgAmp;  % cross-correlation, nt2*nAlphas x nRefs*nVols
        
        mImg=ri.softMask.*img;  % untransformed image
        I2=mImg(:)'*mImg(:);    % image power
        
        R2x=repmat(imgAmp^2*R2(:,iMi)',nTA,1);
        
        diffIR=I2-2*ccIR+R2x;  % squared difference between img and ref, nt2*nAlphas, nRefs*nVols
        
        % Create the WC log priors, which will be nt2*nAlphas x nRefs in
        % size.
        wcPriors1vol=zeros(nTA,nRefs);
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
        wcPriors=zeros(nTA,nRefs,nVols);
        for iVol=1:nVols     % expand width by nVols
            wcPriors(:,:,iVol)=wcPriors1vol+log(ri.pVols(iVol));
        end;
        wcPriors=reshape(wcPriors,nTA,nRV);
        
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
        %             cMeans=zeros(nw^2,nRV,'single');
        %              cMeans=f1(cMeans,eImg_ta_rv*imgAmp,iImg);
        localCMeans(:,:,iMi)=eImg_ta_rv*imgAmp;
        cMeans=cMeans+localCMeans;
        %             ctfos=repmat(ctfo,1,1,nRV);
        %             fclassMeans=fclassMeans+fft2(reshape(cMeans,nw,nw,nRV)).*ctfos*imgAmp;
        %             classNorms=classNorms+kron(ctf(:).^2*imgAmp^2,prv_X(:)'); % nw^2 x nRV
        localCNorms(:,iMi)=imgAmp^2*prv_X(:);  % nMi x nRV
        cNorms=cNorms+localCNorms;
        %             cNorms=f2(cNorms,imgAmp^2*prv_X(:));  % nMi x nRV
        localAmps=ccIR./R2x;
        imgAmps(iImg)=localAmps(:)'*ptarv_X(:)/imgAmp; % divide by imgAmp
        var=diffIR/nPixEff;
        eVar(iImg)=var(:)'*ptarv_X(:);
        
        %             Accumulate the class means
        %             eImg_ta_rvu=eImg_ta_rvu+eImg_ta_rv;
        %             eNorm_rvu=eNorm_rvu+prv_X(:);  % nRefs*nVols x 1
        
        %             pta_Xs(:,iImg)=pta_X;
        pt_Xs(:,iImg)=pt_X;
        
    end; % if
end; % parfor iImg
%     Operate with CTF on classes for the defocus group
%     cMeans=zeros(nw,nw,nRef*nVols,'single');
%     cMeans(pix,:)=eImg_ta_rvu;
%     cMeans=reshape(cMeans,nw,nw,nRV);
%     classMeans=classMeans+real(ifft2(fft2(cMeans).*ctfos));
%     classNorms=classNorms+kron(ctf(:),eNorm_rvu'); %  nw^2 x nRefs*nVols
%     Accumulate the classMeans and classNorms using the ctfs
fclassMeans=zeros(nw,nw,nRV,'single');  % will be complex
classNorms=zeros(nw*nw,nRV,'single');

for iMi=1:nMi
    ctf=si.ctfs(:,:,iMi);
    ctfo=ifftshift(ctf);
    parfor k=1:nRV
        % %             cMean=zeros(nw,nw);
        % %             cMean(pix)=cMeans(:,k,iMi);
        cMean=reshape(cMeans(:,k,iMi),nw,nw);
        fclassMeans(:,:,k)=fclassMeans(:,:,k)+fft2(cMean).*ctfo;
    end;
    classNorms=classNorms+kron(ctf(:).^2,cNorms(:,iMi)'); % nw^2 x nRV
end;
classMeans=reshape(real(ifft2(fclassMeans)),nw,nw,nRefs,nVols);
classNorms=reshape(classNorms,nw,nw,nRefs,nVols);

roi.imgAmps=imgAmps;
pTrans=reshape(pt_Xs,nt,nt,nImgs);
roi.varNs=eVar;
roi.pVols=pVols;
roi.logPX=logPX;
pRefs=reshape(pRefs,nRefs,nVols,nImgs);


% function cM=AccumCMeans(cM,val,i)
%     iM=si.imgIndex(i);
%     cM(:,:,iM)=cM(:,:,iM)+val;
% end
%
% function cM=AccumCNorms(cN,val,i)
%     iM=si.imgIndex(i);
%     cN(:,:,iM)=cN(:,:,iM)+val;
% end

end
