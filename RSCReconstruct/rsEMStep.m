function [moi,classes,norms,probs] = rsEMStep(moi,ri,refs,si,imgs)
%%  Abstraction of the main EM
    
% INPUT:
% imgs
% refs
% -map refs to angles
% -priors on refs
% 
ctfs
ctfIndices(img)
% RunInfo
%     ri.pixA          =   3  % Working angstroms per pixel/voxel in reconstruction, may differ from micrographs
%     ri.n             =  64  % Working box size
%     ri.fWork         =  0.2 % Working bandwidth
%     ri.membraneOffset= -13  % z-coordinate of membrane center in the map.  Negative for AMPAR, positive for BK.
%     ri.maxTranslation=   7 % size of the cross-correlation search is 2*maxTranslation+1
%     ri.alphas
%     ri.refBetas
%     ri.refGammas
%     ri.symmetry      =   2  % We're assuming C symmetry only
%     ri.refIndices    int16angleSteps(2) x angleSteps(3)
% imgAlphaFlags
% imgRefFlags
% imgAmps   

% Model Info
% moi:
%     aI
%     sigmaN
%     sigmaC
%     sigmaG
%     vol
% 
% angles:
%     .mins
%     .steps
%     .num


% output:
% accums
% class means
% class norms
% probs .logP
%       .booleans/pointers for angles
%       .logPs
% 
%       
% p(trans)(img)
% logp(alpha,ref)(img)
% moi: aI
%   sigmaN
%   sigmaC
% 
% probs.logPBeta

% Var naming
% p(t,a|X,r) is named pta_rImg
% order of variables is t,a,r (translation, alpha, reference)

%%  Abstraction of the main EM

nImgs=size(imgs,3);
nAllRefs=size(refs,3);
vesRs=si.vesRs;
ctfs=si.ctfs;
ctfInds=si.ctfInds;

    for iImg=1:nImgs
        alphaInds=find(imgAlphaFlags(:,iImg));
        nAlphas=numel(alphaInds);

        refInds=imgRefFlags(:,iImg);
        nRefs=numel(refInds);
        
        accumImages=zeros(n^2,nRefs);
        accumPars=zeros(npars,nRefs);
        logPr_Img=-inf(nRefs,1);
%         alphas=ri.angleMin(1)+(alphaInds-1)*ri.angleSteps(1);
        alphas=ri.alphas(alphaInds);
        pa_Img=zeros(nAlphas,1); % =pa_img
        par_Img=zeros(nAlphas,nRefs);
        ptar_Img=zeros(n*n,nAlphas,nRefs);
        funrot=fft2(rsRotateImage(imgs(:,:,iImg),alphas));
        
        ctf=ctfs(:,:,ctfIndex(iImg));
        vesR=vesRs(iImg);
        y0=y0s(iImg);
        oldBeta=inf;
            for k=1:nRefs
                refInd=refInds(k);
                beta=ri.betas(refInd);
                if beta~=oldBeta
                    logProbWC=ComputeLogProbWC2(n,moi.sigmaC,moi.sigmaG,vesR,y0,alphas,ri.dAlpha,beta,ri.dBeta);
                    oldBeta=beta;
                end;
                cref=real(ifftn(fftn(refs(:,:,k)).*ctf)).*msk;
                    [logP, accumImg, accumPar, pAlphas, pTrans]...
                        =rsOneImgOneRef6(funrot,logProbWC,imgAmps(iImg),ri.sigmaN,cref);
                    accumImages(:,k)=accumImg(:);                % n^2 x nRefs
                    logPr_Img(k)=logP;
                    accumPars(:,k)=accumPar;
                    pa_rImg(:,k)=pAlphas'; % p(alpha|X,ref)      nAlphas x nRefs
                    pt_arImg(:,:,k)=reshape(pTrans,n*n,nAlphas);  % n^2 x nalpha x nRefs
                end;
            end; % k loop over refs
%             logPr_Img=logPr_Img+log(fmodel);
            maxLogP=max(logPr_Img); %max of k elements
            sclPr=exp(logPr_Img-maxLogP); %k elements
            normPr=sum(sclPr);  % sum_ref: p(X|ref)=p(X,ref) (ref evenly distributed)
            pr_img=sclPr/norm;  % normalized p(k|Img) or p(ref_k|X)=P(X,ref_k)/sum_k(P(X,ref_k)
            par_img=pAlphasR*diag(normPr); % nAlphas*nRefs  x nRefs*nRefs= nAlphas * nRfes
            pa_img=pAlphasR*normPr; %p(alpha, X,ref)=p(alpha|X,ref)*p(X|ref)  nAlphas*nRefs  x nRefs*1  =nAlphas*1
???         
            temp=PAlphas(:,iImg);
            temp(IndexA(:,iImg))=normPa/(sum(normPa)); % nalpha
            PAlphas(:,iImg)=temp; % nAlphas,* nImg
            % to get p(t|X):
            tempt=normPaR/sum(normPa(:));% nalpha*nref (6) P(alpha,ref|X)=P(X,alpha,ref)/sum_alpha,ref(P(X, alpha,ref))
            trans=reshape(pTransR,n*n,nAlphas*nRefs)*tempt(:);% (7) P(t,alpha,ref|X)=P(t|X,alpha,ref)*P(alpha,ref|X)
            % (8) p(t|X)=sum_alpha,ref(P(t,alpha,ref|x))
            nPkIndex=find(normPk>0.0001); % index of references for each image in this iteration
            
            aI(iImg)=accumPars(2,:)*normPk;
            LlI(iImg)=log(norm)+maxLogP; % per image
            accumSums=accumSums+accumPars*normPk;
            accumsI(:,iImg)=accumPars*normPk;
            pixNormPk=repmat(normPk',n^2,1);

% %             ctAccumImg=real(ifftn(fftn(accumImg).*ctf));
% %                     accumImages(:,k)=ctAccumImg(:);

            classSums1=classSums1+accumImages.*repmat(normPk',n^2,1)*mod(iImg,2);  % npixels x nrefs
            classNorms1=classNorms1+ctf(:).^2*normPk'*mod(iImg,2);
            classSums2=classSums2+accumImages.*repmat(normPk',n^2,1)*(~mod(iImg,2));
            classNorms2=classNorms2+ctf(:).^2*normPk'*(~mod(iImg,2));
    end; % iImg loop over images
    toc
