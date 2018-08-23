
    
    modelInfo
    mo.noiseSD    % noise standard deviation
    mo.clickSD    % in pixels
    mo.bobMean
    mo.pVols(k)      % probability of the various 3D volumes
    mo.likelihood
    
    
    input:
    imgs
    
    refs
    -map refs to angles
    -priors on refs
    ctfs
    ctf(img) indices
    a(img) vals
    logp(alpha,ref)(img)
    pars:
    sigmaN
    sigmaC
    activeListThreshold
    
    output:
    accums
    class means
    class norms
    logPs
    p(trans)(img)
    logp(alpha,ref)(img)
    a(img)
    
    probs.logPBeta
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%% new code %%%%%%%%%%%%%%%
    
    % realMask=...
    
    % Common to all
    rfi.refs=refs;
    % rfi.angles=angs; (nProjAngles x 3); nProjAngles = nBetas x nGammas
    % rfi.alphas=...; (nalphas x 1)
    % rfi.betas=...;
    rfi.priors=refPriors;
    % rfi.dAlpha=...
    % rfi.dBeta=...
    
    % logPs.angleMap: uint8 nalphas x nbetas x ngammas

    
    
    
    for iImg=1:numImgs
        imi=imis(iImg);
        % For this image
        if aI(iImg)>0  % don't use it if its amplitude is zero.
            
            if iIter==1
                imi.amp=aI(iImg);
                imi.vesR=si.rVesicle(iImg);
                imi.yClick=si.yClick(iImg);
                imi.logProbWC=zeros(nAlphas,nBetas,'single');
                
                for k=1:rfi.nBetas % calculate p(w,c) prior
                    beta=betas(k);
                    imi.logProbWC(:,k)=ComputeLogProbWC2(1,pars.sigmaC,pars.sigmaG,...
                        si.rVesicle(iImg),si.yClick(iImg),alphas,dAlpha,beta,dBeta);
                end;
                
                %     Assign the angleMap (nAlphas x nBetas x nGammas) to the priors
                logPs.angleMap=repmat(uint8(floor(logProbWC/logPs.alphaBetaScale+255)),[1 1 nGammas]);
                imi.funrot=zeros([n 0],'single');
                imi.unrotAlphaInds=zeros(1,0,'int16');
            end;

            % Assign the activeList based on the priors
%             This can 
            activeList=logPs(iImg).angleMap>=pars.activeListThreshold; % nalphas x nbetas x ngammas
            % Figure out which alphas we have to evaluate
            alphaIndList=find(any(reshape(activeList,nAlphas,nProjAngles),2));
%             Downsample and filter the image
%             
            [filtImg,alphaInds,funrotImgs]=GetUnrotImgs(filtImg,fmasko,rfi.alphas,alphaIndList,alphaInds,funrotImgs);

            function imi=FilterAndGetUnrotImgs(imi,fmask,alphas,alphaIndList);
% Expand the funrotImgs array to include all the entries in the
% alphaIndList. Put the filtered image and a complete set of unrot images
% into the imi.funrot and imi.unrotAlphaInds fields.

            
            
            
            nAlphasToEval=numel(alphaIndList);
            if nAlphasToEval<1
                error('No alphas to evaluate');
            end;
            addEvalList=[];
            for i=1:nAlphasToEval
                if ~any(imi.alphaIndList==alphaIndList(i))
                    addEvalList=[addEvalList allphaIndList(i)];
                end;
            end;
            imi.funrot(:,:,end+1:end+numel(addEvalList))=UnrotateImages(img,rfi.alphas(addEvalList));
            imi.unrotAlphaInds=[unrotAlphaInds addEvalList];
            
            % Fill in the rest of the image info structure
            imi.img=imgs(:,:,iImg);
            imi.vesR=si.rVesicle(iImg);
            imi.y0=si.y0s(iImg);
            imi.ctf=imi.amp*ctfs(:,:,ctfIndex(iImg));  % weight by amplitude here.
% structure for rotated images
            iri.funrot=single(0);
            iri.unrotAlphaInds=int16(1,0);
            
%%%%%%%%%%%% Abstraction of EM for one image
            [logP, iri, accums, probs]=EMOneImage(imi, iri, rfi, activeList, refIndices, pars);
            % imgInfo.funrotImgs
            % imgInfo.alphas
            % imgInfo.ctf
            % imgInfo.realMask
            % imgInfo.a
            % refInfo.refs
            % refInfo.angles
            % refInfo.priors
            % activeList = logical(nalphas,nrefs)
%             probs are returned only for some references, indexed by
%             probs.refIndices.  However, all alpha values are evaluated.
            nPars=5;
            n=size(imi.img);
            npix=prod(n);
            [nAllAlphas nAllBetas nAllGammas]=size(activeList);
            activeRefs=any(activeList);
            probs.RefIndices=int16(sum(activeRefs),1);
%             refIndices=find(activeRefs);
%             nRefs=sum(activeRefs);
            betaInds=find(any(activeRefs,3));
            nBetas=numel(betaInds);
            
            accumImages=zeros(npix,nRefs,'single');
            logPs=zeros(nRefs,1);
            accums=zeros(nPars,nRefs);
            
            pAlphasR=zeros(nAllAlphas,nRefs);
            pTransR=zeros(npix,nAllAlphas,nRefs);
            
            iRef=0;
            
            for iBeta=1:nBetas
                %     Find all the alphas active at any of the gamma values, and compute
                %     the prior
                betaInd=betaInds(iBeta);
                beta=rfi.betas(betaInd);
                activeAlphas=find(any(activeList(:,betaInd,:),3));
%                 Prior from sigmaC and sigmaG
                logProbWC=ComputeLogProbWC2(n,pars.sigmaC,pars.sigmaG,imi.vesR,imi.y0,...
                    rfi.alphas(activeAlphas),rfi.dAlpha,beta,rfi.dBeta);
                %     Evaluate the active gammas, to enumerate all the refr
                activeGammas=find(any(activeList(:,betaInds(iBeta),:)));
                nGammas=numel(activeGammas);
                for iGamma=1:nGammas
                    gammaInd=activeGammas(iGamma);
                    iRef=iRef+1;
                    refInd=sub2ind([nAllBetas nAllGammas],betaInd,gammaInd);  % reference index
                    probs.refIndices(iRef)=refInd;
                    alphaInds=find(activeList(:,refInd));
                    %     alphas=imi.alphas(alphaInds);
                    unrots=imi.funrot(:,:,alphaInds);
                    beta=rfi.angs(iRef,2);
                    
                    cref=real(ifftn(fftn(refs(:,:,iRef)).*imi.ctf)).*realMask;
%                     Main computation
                    [logP, fAccumImg, accumPar, pAlphas, pTrans]...
                        =rsOneImgOneRef6(funrot,logProbWC,imi.amp,pars.sigmaN,cref);
%                     Filter the accumulated image by the ctf.
                    accumImgC=real(ifftn(ifftshift(imi.ctf).*fAccumImg));
                    accumImages(:,iRef)=accumImgC(:);  %  E{img}_alpha,trans | ref
                    logPs(iRef)=logP;                 %  p(ref | X)
                    accumPars(:,iRef)=accumPar;          %  E{sums}_alpha,trans
                    pAlphasR(alphaInds,iRef)=pAlphas';%  p(alpha|X,ref)      nAlphas x nRefs
                    pTransR(:,iRef)=pTrans*pAlphasR;  %  p(trans|ref);
                    pTransAR(:,:,iRef)=pTrans$$$$$$$$$$$
                    %           p(trans | alpha, ref)
                end; % loop over gammas
            end; % k loop over betas

%             get p(ref_k|X)
            logPs=logPs+log(rfi.priors(refIndices'));
            maxLogP=max(logPs); %max of nRefs elements
            sclPk=exp(logPs-maxLogP); % nRefs elements
            norm=sum(sclPk);  % sum_ref: p(X|ref)=p(X,ref) (ref evenly distributed)
            normPk=sclPk/norm;  % column vector, p(ref_k|X)=P(X,ref_k)/sum_k(P(X,ref_k)

%             get p(alpha,ref|X) and p(alpha | X)
            normPaR=pAlphasR*diag(normPk); % joint p(alpha,ref|X) = p(alpha|ref,X)*p(ref|X).  nAlphas*nRefs
            normPa=sum(normPaR,2);         % p(alpha)
%             normPa=pAlphasR*normPk; %p(alpha| X,ref)=p(alpha|X,ref)*p(ref|X)  nAlphas*nRefs  x nRefs*1  =nAlphas*1
            %             =sum(normPaR,2)
            accums.sums=accumSums+accumPars*normPk;
            
            normPkImage=repmat(normPk',npix,1);
%             accums.classSums=accumImages.*normPkImage;  % npix x nRefs
            accums.classSums=accumImages*diag(normPk);
            
            normImage=real(ifftn(ifftshift(imi.ctf)));
            accums.classNorms(:,iImg)=normImage(:)*normPk';  %npix x nRefs

            probs.logPRefs=log(normPk);
            probs.logPAlpha=log(normPa);
            probs.logPTrans=
            
% We wish to create the following outputs:
% logP  The likelihood for the one image
% accums.pars
% accums.classSums            
% accums.classNorms
% probs.logPTrans   ln p( trans | X )
% probs.logPAlpha   ln p( alpha | X )
% probs.logPRefs    ln p( refs  | X )
% probs.logPMap     ln p(alpha beta gamma) in bytes


            
            
            
            
            
            for iImg=1:nImgs
                if (ctfIndex(iImg))
                    accumImages=zeros(n^2,nRefs);
                    accumPars=zeros(npars,nRefs);
                    logPs=-inf(nRefs,1);
                    alphas=alphas0(IndexA(:,iImg));
                    nAlphas=numel(alphas);
                    accumPAlphas=zeros(nAlphas,1);
                    %PAlphasB=zeros(nAlphas, nRefBetas,nRefGammas);
                    PAlphasR=zeros(nAlphas,nRefs);
                    pTransR=zeros(n*n,nAlphas,nRefs);
                    funrot=fft2(rsRotateImage(imgs(:,:,iImg),alphas));
                    
                    ctf=ctfs(:,:,ctfIndex(iImg));
                    vesR=vesRs(iImg);
                    
                    nPkIndex=nPkIndice(:,iImg);
                    %     Loop over references
                    for k=1:nRefs
                        if find(nPkIndex==k)%k<=nPkIndex(k)
                            if mod(k,nRefs0)
                                beta=angs(mod(k,nRefs0),2);
                            else
                                beta=angs(nRefs0,2);
                            end
                            logProbWC=ComputeLogProbWC2(n,sigmaC,sigmaG,vesR,y0s(iImg),alphas,dAlpha,beta,dBeta);
                            %logProbWC=ComputeLogProbWC2(n,sigmaC,sigmaG,vesR+aS(4),y0s(iImg),alphas,dAlpha,beta,dBeta);
                            %         logProbWC=zeros(n,n,numel(alphas));
                            %        cref=real(ifftn(fftn(refs(:,:,k)).*ctfs(:,:,iImg)));
                            cref=real(ifftn(fftn(refs(:,:,k)).*ctf)).*msk;
                            [logP accumImg accumPar pAlphas pTrans]=rsOneImgOneRef5(funrot,logProbWC,aI(iImg),sigmaN,cref);
                            %loP: one element; pALphas:nAlphas pTrans: n*ns
                            ctAccumImg=real(ifftn(fftn(accumImg).*ctf));
                            accumImages(:,k)=ctAccumImg(:);
                            logPs(k)=logP;
                            accumPars(:,k)=accumPar;
                            if sum(isnan(accumPar))
                                disp([num2str([iImg,k,accumPar(:)]),' accumPar']);
                            end;
                            %accumPAlphas=accumPAlphas+pAlphas;
                            PAlphasR(:,k)=pAlphas'; % p(alpha|X,ref)      nAlphas * nRefs
                            %PAlphasB(:,fix(k,nRefGammas),mod(k,nRefGammas)+1)=pAlphas'; %  each ref and alpha
                            pTransR(:,:,k)=reshape(pTrans,n*n,nAlphas);  % n*n nalpha nref
                        end;
                    end; % k loop over refs
                    logPs=logPs+log(fmodel);
                    maxLogP=max(logPs); %max of k elements
                    sclPk=exp(logPs-maxLogP); %k elements
                    % %     [norm1,norm2]=sum(reshape(sclPk,nRefs0,nmodel),2);
                    norm1=sum(sclPk(1:nRefs0));
                    norm2=sum(sclPk(1+nRefs0:nRefs));
                    norm=sum(sclPk);  % sum_ref: p(X|ref)=p(X,ref) (ref evenly distributed)
                    normPk=sclPk/norm;  % normalized p(k|Img) or p(ref_k|X)=P(X,ref_k)/sum_k(P(X,ref_k)
                    normPaR=PAlphasR*diag(normPk); % nAlphas*nRefs  x nRefs*nRefs= nAlphas * nRfes
                    normPa=PAlphasR*normPk; %p(alpha, X,ref)=p(alpha|X,ref)*p(X|ref)  nAlphas*nRefs  x nRefs*1  =nAlphas*1
                    %(4) P(X,alpha)=sum_ref(P(X,alpha,ref))     nAlphas
                    %norma=sum(normPa,2);
                    %accumPAlphasRI(:,:,iImg)=normPa./repmat(norma, [1,nRefs]);
                    temp=PAlphas(:,iImg);
                    temp(IndexA(:,iImg))=normPa/(sum(normPa)); % nalpha
                    PAlphas(:,iImg)=temp; % nAlphas,* nImg
                    %PAlphas(IndexA(:,iImg),iImg)=norma/sum(norma); % nAlphas,* nImg
                    %P(alpha|X)=P(X,alpha)/sum_alpha(P(X,alpha))
                    
                    %accumPalphasBI(:,:,iImg)=sum(normPa())
                    %     for iBeta=1:nRefBetas
                    %         accumBeta(iBeta,iGamma)=sum(normPk((0+iBeta):nRefGammas:(nRef-nRefGammas+iBeta)));
                    %     end;
                    
                    % to get p(t|X):
                    tempt=normPaR/sum(normPa(:));% nalpha*nref (6) P(alpha,ref|X)=P(X,alpha,ref)/sum_alpha,ref(P(X, alpha,ref))
                    trans=reshape(pTransR,n*n,nAlphas*nRefs)*tempt(:);% (7) P(t,alpha,ref|X)=P(t|X,alpha,ref)*P(alpha,ref|X)
                    % (8) p(t|X)=sum_alpha,ref(P(t,alpha,ref|x))
                    %        if iImg<71
                    %            imacs(reshape(trans,n,n));
                    %             colorbar;
                    %             title(num2str(iImg));
                    %             drawnow;
                    %
                    %             frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                    %             writeVideo(writerObj, frame);
                    %
                    %        end;
                    nPkIndex=find(normPk>0.0001); % index of references for each image in this iteration
                    
                    if iiter==1
                        nPkIndex0=zeros(nRefs,1);
                        nPkIndex0(1:numel(nPkIndex))=nPkIndex;
                        %nPkIndice(:,iImg)=0;
                        nPkIndice(:,iImg)=nPkIndex0;
                        
                        %if (max(PAlphas)./(min(PAlphas)+1E-10))>0.001
                    end;
                    aI(iImg)=accumPars(2,:)*normPk;
                    %accumL=norm*exp(maxLogP)+accumL;
                    LlI1(iImg)=log(norm1)+maxLogP;
                    LlI2(iImg)=log(norm2)+maxLogP;
                    
                    LlI(iImg)=log(norm)+maxLogP; % per image
                    accumSums=accumSums+accumPars*normPk;
                    accumsI(:,iImg)=accumPars*normPk;
                    pixNormPk=repmat(normPk',n^2,1);
                    
                    classSums1=classSums1+accumImages.*repmat(normPk',n^2,1)*mod(iImg,2);  % npixels x nrefs
                    classNorms1=classNorms1+ctf(:).^2*normPk'*mod(iImg,2);
                    classSums2=classSums2+accumImages.*repmat(normPk',n^2,1)*(~mod(iImg,2));
                    classNorms2=classNorms2+ctf(:).^2*normPk'*(~mod(iImg,2));
                    % %     Show the best aligned image
                    %     subplot(221);
                    %     [mxVal mxK]=max(normPk);  % show the reference with the best alpha
                    %     imacs(reshape(accumImages(:,mxK),n,n));
                    %     title(num2str([iImg mxK]));
                    %
                    %     subplot(222);
                    %     plot(alphas,accumPAlphas/iImg);
                    %     ylabel('pAlphas');
                    %     xlabel('Alpha');
                    %     drawnow;
                end;  % if ctfIndex
                %     if iImg==70
                %         close(writerObj); % Saves the movie.
                %     end;
            end; % iImg loop over images
            toc
            %%  Display the results
            
            if iiter==1
                IndexA=PAlphas>0.001;
            end;
            maxlogPI=max(LlI);
            accumL=sum(exp(LlI-maxlogPI)); % all imges
            accumL1=sum(exp(LlI1-maxlogPI));
            accumL2=sum(exp(LlI2-maxlogPI));
            fmodelI=[accumL1,accumL2]/(accumL1+accumL2)
            fmodel=reshape(repmat(fmodelI,nRefs0,nmodel),nRefs,1);
            accumLl=log(accumL)+maxlogPI
            figure(1);
            %
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            %%  Abstraction of the main EM
            [moi,classes,norms,probs] = EMStep(moi,si,img,ctfs,refs,angles)
            
            nImgs=size(imgs,3);
            nRefs=size(refs,3);
            
            
            for iImg=1:nImgs
                accumImages=zeros(n^2,nRefs);
                accumPars=zeros(npars,nRefs);
                logPs=-inf(nRefs,1);
                alphas=alphas0(IndexA(:,iImg));
                nAlphas=numel(alphas);
                accumPAlphas=zeros(nAlphas,1);
                %PAlphasB=zeros(nAlphas, nRefBetas,nRefGammas);
                pAlphasR=zeros(nAlphas,nRefs);
                pTransR=zeros(n*n,nAlphas,nRefs);
                funrot=fft2(rsRotateImage(imgs(:,:,iImg),alphas));
                
                ctf=ctfs(:,:,ctfIndex(iImg));
                vesR=vesRs(iImg);
                
                operating on 1 image:
                
                refSelection = logical(nbetas, nalphas)
                
                refinfo:
                [angs
                    refs
                    priors on refs]
                
                imginfo
                [unrotImgs
                    ctf
                    vesR
                    a]
                pars
                (img and refs have been pre-masked; classes have to be masked)
                
                outputs
                
                accums
                a
                
                logP
                
                classes:
                [class means
                    class norms]
                
                imageProbs:
                [p(trans)(img)
                    logp(alpha,ref)(img)]
                
                
                
                
                
                nPkIndex=nPkIndice(:,iImg);
                %         for iModel=1:nModels
                %     Loop over references
                for k=1:nRefs
                    if find(nPkIndex==k)%k<=nPkIndex(k)
                        if mod(k,nRefs0)
                            beta=angs(mod(k,nRefs0),2);
                        else
                            beta=angs(nRefs0,2);
                        end
                        logProbWC=ComputeLogProbWC2(n,sigmaC,sigmaG,vesR,y0s(iImg),alphas,dAlpha,beta,dBeta);
                        cref=real(ifftn(fftn(refs(:,:,k)).*ctf)).*msk;
                        [logP, accumImg, accumPar, pAlphas, pTrans]=rsOneImgOneRef6(funrot,logProbWC,aI(iImg),sigmaN,cref);
                        %loP: one element; pALphas:nAlphas pTrans: n*ns
                        accumImages(:,k)=accumImg(:);                % n^2 x nRefs
                        logPs(k)=logP;
                        accumPars(:,k)=accumPar;
                        pAlphasR(:,k)=pAlphas'; % p(alpha|X,ref)      nAlphas x nRefs
                        pTransR(:,:,k)=reshape(pTrans,n*n,nAlphas);  % n^2 x nalpha x nRefs
                    end;
                end; % k loop over refs
                logPs=logPs+log(fmodel);
                maxLogP=max(logPs); %max of k elements
                sclPk=exp(logPs-maxLogP); %k elements
                % %     [norm1,norm2]=sum(reshape(sclPk,nRefs0,nmodel),2);
                norm1=sum(sclPk(1:nRefs0));
                norm2=sum(sclPk(1+nRefs0:nRefs));
                norm=sum(sclPk);  % sum_ref: p(X|ref)=p(X,ref) (ref evenly distributed)
                normPk=sclPk/norm;  % normalized p(k|Img) or p(ref_k|X)=P(X,ref_k)/sum_k(P(X,ref_k)
                normPaR=pAlphasR*diag(normPk); % nAlphas*nRefs  x nRefs*nRefs= nAlphas * nRfes
                normPa=pAlphasR*normPk; %p(alpha, X,ref)=p(alpha|X,ref)*p(X|ref)  nAlphas*nRefs  x nRefs*1  =nAlphas*1
                
                temp=PAlphas(:,iImg);
                temp(IndexA(:,iImg))=normPa/(sum(normPa)); % nalpha
                PAlphas(:,iImg)=temp; % nAlphas,* nImg
                %PAlphas(IndexA(:,iImg),iImg)=norma/sum(norma); % nAlphas,* nImg
                %P(alpha|X)=P(X,alpha)/sum_alpha(P(X,alpha))
                
                %accumPalphasBI(:,:,iImg)=sum(normPa())
                %     for iBeta=1:nRefBetas
                %         accumBeta(iBeta,iGamma)=sum(normPk((0+iBeta):nRefGammas:(nRef-nRefGammas+iBeta)));
                %     end;
                
                % to get p(t|X):
                tempt=normPaR/sum(normPa(:));% nalpha*nref (6) P(alpha,ref|X)=P(X,alpha,ref)/sum_alpha,ref(P(X, alpha,ref))
                trans=reshape(pTransR,n*n,nAlphas*nRefs)*tempt(:);% (7) P(t,alpha,ref|X)=P(t|X,alpha,ref)*P(alpha,ref|X)
                % (8) p(t|X)=sum_alpha,ref(P(t,alpha,ref|x))
                %        if iImg<71
                %            imacs(reshape(trans,n,n));
                %             colorbar;
                %             title(num2str(iImg));
                %             drawnow;
                %
                %             frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                %             writeVideo(writerObj, frame);
                %
                %        end;
                nPkIndex=find(normPk>0.0001); % index of references for each image in this iteration
                
                if iiter==1
                    nPkIndex0=zeros(nRefs,1);
                    nPkIndex0(1:numel(nPkIndex))=nPkIndex;
                    %nPkIndice(:,iImg)=0;
                    nPkIndice(:,iImg)=nPkIndex0;
                    
                    %if (max(PAlphas)./(min(PAlphas)+1E-10))>0.001
                end;
                aI(iImg)=accumPars(2,:)*normPk;
                %accumL=norm*exp(maxLogP)+accumL;
                LlI1(iImg)=log(norm1)+maxLogP;
                LlI2(iImg)=log(norm2)+maxLogP;
                
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
                % %     Show the best aligned image
                %     subplot(221);
                %     [mxVal mxK]=max(normPk);  % show the reference with the best alpha
                %     imacs(reshape(accumImages(:,mxK),n,n));
                %     title(num2str([iImg mxK]));
                %
                %     subplot(222);
                %     plot(alphas,accumPAlphas/iImg);
                %     ylabel('pAlphas');
                %     xlabel('Alpha');
                %     drawnow;
            end;  % if ctfIndex
            %     if iImg==70
            %         close(writerObj); % Saves the movie.
            %     end;
        end; % iImg loop over images
        toc
        %%  Display the results
        
        if iiter==1
            IndexA=PAlphas>0.001;
        end;
        maxlogPI=max(LlI);
        accumL=sum(exp(LlI-maxlogPI)); % all imges
        accumL1=sum(exp(LlI1-maxlogPI));
        accumL2=sum(exp(LlI2-maxlogPI));
        fmodelI=[accumL1,accumL2]/(accumL1+accumL2)
        fmodel=reshape(repmat(fmodelI,nRefs0,nmodel),nRefs,1);
        accumLl=log(accumL)+maxlogPI
        figure(1);
        %disp(accumSums/nImgs);  % accumulators
        subplot(222);
        plot(accumsI(2,:)); title('Amplitudes');
        xlabel('Image no.');
        set(gca,'fontsize',14);
        
        figure(2);
        SetGrayscale;
        clf;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%need editting for
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fsc real
        classSums=classSums1+classSums2;
        classNorms=classNorms1+classNorms2;
        clsSum=reshape(classSums,n,n,nRefs);
        clsNorm=reshape(classNorms,n,n,nRefs);
        k=.2;
        classMeans=real(ifft2(fft2(clsSum)./(k+clsNorm)));
        ixs=.5:nRefGammas+.5;
        iys=.5:nRefBetas+.5;
        imacs(ixs,iys,ImageArray(classMeans,0,n,1,nRefGammas,nmodel*nRefBetas));
        set(gca,'xtick',1:nRefGammas,'XTickLabel',num2str(refGammas'));
        xlabel('Gamma');
        set(gca,'ytick',1:nRefBetas,'YTickLabel',num2str(refBetas'));
        ylabel('Beta');
        
        figure(1);
        subplot(223);
        semilogy(sum(classNorms),'-');
        axis([0 nRefs 1e-8 inf]);
        set(gca,'fontsize',14);
        % hold on; % Mark the references of the original images
        % semilogy(refPtrs,classNorms(refPtrs),'r+');
        % hold off;
        xlabel('Reference no.');
        ylabel('Class norm');
        
        subplot(224);
        [mxVal mxK]=max(sum(classNorms));
        bestImg=reshape(classSums(:,mxK),n,n);
        imacs(bestImg);
        title('Strongest class mean');
        set(gca,'fontsize',14);
        aS=accumSums/sum(logical(find(ctfIndex>0)));
        %aS=accumSums/nImgs;
        aS(1)=aS(1)*n^2/(msk(:)'*msk(:));
        % aS(3)=aS(3)-aS(4)^2;
        %theoret=[sigmaN^2 a sigmaC^2 b0 0];
        if ~Fake
            theor_a=aS(2)*si.sVesicle/mean(si.sVesicle);
            ctfIndex(logical((aI./theor_a>1.8)+(aI./theor_a<0.55)))=0;
        else
            ctfIndex(logical((aI<1/10*aS(2))+(aI>10*aS(2))))=0; % neglect imgs with absurd a
        end;
        textCells={'noise var  '
            'a          '
            'click var  '
            'b0         '
            '           '};
        for i=1:numel(accumSums)
            disp([textCells{i} num2str([aS(i) theoret(i)])]);
        end;
        %%
        % clsNorm=repmat(ifftshift(fuzzymask(n,2,0.4*n,.05*n)),[1 1 nRefs]);
        % rclsNorm=zeros(n,n,nRefs);
        % clsSum=refs;
        rclsNorm=zeros(n,n,nRefs);
        for i=1:nRefs
            rclsNorm(:,:,i)=fftshift(ifft2(clsNorm(:,:,i)));
        end;
        symmetry=2;
        %%
        %[reconVol normVol]=rsDoReconstructionSymmetry(clsSum,rclsNorm,angs,symmetry);
        for imodel=1:nmodel
            figure(2+(iiter-1)*nmodel+imodel);
            % reconstruction
            k=1;
            %mapUM=rsNormalizeReconstruction(reconVol,normVol,k);
            
            
            % FSC
            %fsc=FSCbyhalves(classSums1,classSums2,classNorms1,classNorms2,n,nRefs,angs);
            %fsc=FSCorr(vol,mapx);  % Compare with crystal structure
            
            range=nRefs0*(imodel-1)+1:nRefs0*imodel;
            [fsc,mapUM]=VFSCbyhalvesM(k,symmetry,classSums1(:,range),classSums2(:,range),classNorms1(:,range),classNorms2(:,range),n,nRefs0,angs,refMsks);
            map=mapUM.*Volmsk;
            vol=map*std(vol(:))/std(map(:));  % scale
            ShowSections(vol);
            fscs=[fscs,fsc];
            subplot(339);
            %plot(fscs);
            plot(fscs(:,1+imodel:nmodel:nfsc));
            set(gca,'xtick',1:3:n/2,'XTickLabel',round(n*pixA./(1:3:n/2)'));
            xlabel('resolution (Å)');
            line(pixA*n./[30,30],[0,1]);
            line(pixA*n./[25,25],[0,1]);
            line(pixA*n./[20,20],[0,1]);
            line([0,n/2],[0.5,0.5]);
            ML.V{imodel}(:,:,:,iiter)=vol;
        end;
        %% store the values
        
        ML.sigmaN2(iiter)=aS(1);
        ML.a(iiter)=aS(2);
        ML.sigmaC2(iiter)=aS(3);
        ML.b0(iiter)=aS(4);
        Ml.L(iiter)=accumLl;
    end;
    disp('Done.');
    save('ResultsOldData10Inter.mat','ML');
    
    
