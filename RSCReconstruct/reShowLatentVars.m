function reShowLatentVars(imgs,refs,ri,moi,roi,iter,figs,imIndices)
nImgs=size(imgs,3);
if nargin<8
    imIndices=round([1/nImgs .25 .5 .5+1/nImgs 1]*nImgs);
end;

maxShift=(ri.nTrans-1)/2;
nVols=size(roi.pVols,1);
nRefs=size(refs,3);
alphas=ri.alphas(:,1);
nAlphas=numel(alphas);
nIso=size(ri.alphas,2);
    nt=sqrt(size(roi.pTrans,1));

%         q2=reshape(roi.pRefs,ri.angleN(3),ri.angleN(2),nVols,nImgs);
    betas=ri.angleMin(2)+(0:ri.angleN(2)-1)'*ri.angleStep(2);
    gammas=ri.angleMin(3)+(0:ri.angleN(3)-1)'*ri.angleStep(3);

%%
if figs(1)>0
    figure(figs(1));  % Show the overall latent probability distributions
    SetGrayscale;
    subplot(2,2,1);
    plot(roi.imgAmps);
    xlabel('Image');
    ylabel('Image amplitude');
    title(['Iteration ' num2str(iter)]);
    
%     subplot(2,2,2);
    subplot(4,4,3);
    ixs=-maxShift:maxShift;
    
    pt=reshape(mean(roi.pTrans,2),nt,nt);
    imacs(ixs,ixs,pt.^.2);
    xlabel('X shift');
    ylabel('Y shift');
    
%     subplot(2,2,2);
%     na=size(roi.pAlphas,1)*size(roi.pAlphas,2);
%     if numel(alphas)<na
%         alphas=1:na;
%     end;
    subplot(4,4,7);
    bar(alphas,mean(roi.pAlphas,3));
    xlabel('AlphaI');
    ylabel('P(alpha)');
    
    subplot(2,2,3);
    q=reshape(roi.pRefs,nRefs*nVols,nImgs);
%     imacs(log(max(1e-2,q')));
    imacs(q'.^.2);
    xlabel('Image');
    ylabel('Reference');
    
    ctr=floor([ri.angleN(3) ri.angleN(2)]/2)+1;
    avgP=zeros(ri.angleN(3),ri.angleN(2));
    volP=zeros(nVols,1);
    ix=volP;
    iy=volP;
    for i=1:nImgs
        for j=1:nVols
            [volP(j),ix(j),iy(j)]=max2d(roi.pRefs(:,:,j,i));
        end;
        [mxVal,j0]=max(volP);
        avgP=avgP+circshift(roi.pRefs(:,:,j0,i),ctr-[ix(j) iy(j)]);
    end;

    avgP=avgP/max(avgP(:));
    subplot(224);
    imacs(gammas,betas,avgP.^.2);
    xlabel('Gamma');
    ylabel('Beta');
%     colormap jet
%     colorbar;
end;

% ----------------------------------

    figure(figs(2));  % Show exemplary particles and matches
    imIndices=min(imIndices,nImgs);
    [~,iRefs]=max(reshape(roi.pRefs,nRefs*nVols,nImgs),[],1);
    [~,iVolVs]=ind2sub([nRefs nVols],iRefs);
    nc=5;
    nr=numel(imIndices);  % number of rows
    oColor=[0 .2 .4];  % mask overlay color
    tColor=[.4 .2 0];
    lowThresh=1e-3;
    
    for i=1:nr
        ind=imIndices(i);
        
% First column: show the soft-masked image
        mysubplot(nr,nc,1+nc*(i-1));
        imags(imgs(:,:,ind).*ri.softMask);
        axis equal off;
        if isfield(roi,'inBoundImgs') && ~roi.inBoundImgs(ind)
            inBoundText=' X';
        else
            inBoundText=' ';
        end;
        title(['Image ' num2str(ind) inBoundText]);

% 2.  Show the best-matching reference
        mysubplot(nr,nc,2+nc*(i-1));
        imags(refs(:,:,iRefs(ind)));
        if i==1
                title('Best match');
        end;
%         axis off;
        
% 3.  Create a grayscale image with mask overlay ot the trans probs
        mysubplot(nr,nc,3+nc*(i-1));
        pt=reshape(roi.pTrans(:,ind),nt,nt)';
        pm=reshape(~moi.activeTrans(:,ind),nt,nt)';
        pmsk=pm(:)*oColor; % nt2 x 3
        pmsk=reshape(pmsk,nt,nt,3);
        ptImg=pmsk+repmat(pt/max(pt(:)),1,1,3);  % make rgb grayscale
%         image(ptImg'); axis xy
%         image(-maxShift:maxShift,-maxShift:maxShift,ptImg);
        image(ptImg);
        axis xy off
        if i==1
                title(['Iteration ' num2str(iter)]);
        end;
%         xlabel('X shift');
%         ylabel('Y shift');

% 4.  Plot p(alpha), with rso and iso in different colors.
        mysubplot(nr,nc,4+nc*(i-1));
        pAlpha=roi.pAlphas(:,:,ind);
        msk=moi.activeAlphas(:,:,ind);
        amsk=single(msk);
        amsk(~msk)=NaN;
        plot(alphas,pAlpha);
        hold on;
        plot(alphas,pAlpha.*amsk,'k.','markersize',20);
        hold off;
%         xlabel('Alpha');

% 5.  Plot p(beta,gamma)
        mysubplot(nr,nc,5+nc*(i-1));
        q3=roi.pRefs(:,:,iVolVs(ind),ind)';
        [nBeta,nGamma]=size(q3);
        q3=q3/(max(q3(:)));
        aRVs=moi.activeRVs(:,:,iVolVs(ind),ind)';
        q3Low=~aRVs(:)*tColor;
        q3Img=repmat(q3,1,1,3)+reshape(q3Low,nBeta,nGamma,3);
        image(gammas,betas,q3Img);
%         axis xy
%         xlabel('Gamma');
%         ylabel('Beta');
        
    end;

