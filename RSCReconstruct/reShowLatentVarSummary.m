% reShowLatentVarSummary(ri,roi)
nargin=0;
if nargin<2
    disp('Loading the roi file');
    [name, pathName]=uigetfile('*roi.mat','Select roi file');
    if isnumeric(pathName) % File selection was cancelled
        return
    end;
    cd(pathName);
load(name);  % get the roi file
load('ri.mat'); % get the latest ri
end
%%
ri0=ri;
[nGammas,nBetas,nVols,nIm]=size(roi.pRefs);
if size(ri.angles,1)==nBetas*nGammas
    ri=ri0;
else
    ok=false;
    for iseq=1:size(ri0.nSequence,1)
        n=ri0.nSequence(iseq,1);
        ri=reMakeRunInfoScaled(ri0,n);
        nCurrent=ri.nCurrent;
        
        if size(ri.angles,1)==nBetas*nGammas
            ok=true;
            break;
        end;
    end;
    if ~ok
        error(['No match found for numAngles=' num2str(size(ri0.angles,1)) ' and roi numRefs= ' num2str(nBetas*nGammas)]);
    end;
end;
n=ri.nCurrent;

maxShift=(ri.nTrans-1)/2;
nVols=size(roi.pVols,1);
alphas=ri.alphas(:,1);
nAlphas=numel(alphas);
nIso=size(ri.alphas,2);
nt=sqrt(size(roi.pTrans,1));
nImgs=numel(roi.imgAmps);
nRefs=nBetas*nGammas*nVols;
    betas=ri.angleMin(2)+(0:ri.angleN(2)-1)'*ri.angleStep(2);
    gammas=ri.angleMin(3)+(0:ri.angleN(3)-1)'*ri.angleStep(3);

    subplot(2,3,1);
    plot(roi.imgAmps);
    xlabel('Image');
    ylabel('Image amplitude');
    
    ixs=-maxShift:maxShift;
    
    pt=reshape(mean(roi.pTrans,2),nt,nt);
    imags(ixs,ixs,pt.^.2);
    xlabel('X shift');
    ylabel('Y shift');
    
    subplot(232);
    bar(alphas,mean(roi.pAlphas,3));
    xlabel('AlphaI');
    ylabel('P(alpha)');
    
    subplot(2,3,3);
    nRV=nRefs*nVols;
    q=reshape(roi.pRefs,nRV,nImgs);
    ds=ceil(max(nRV,nImgs)/512);
    imags(max(0,GaussFilt(q',.4/ds)).^.2);
    xlabel('Image');
    ylabel('Reference');
    %%
    subplot(234);
    betaRefs=reshape(sum(roi.pRefs,1),nBetas*nVols,nImgs);
    betaVals=repmat(betas,nVols,1);
    imags(1:nImgs,betaVals,betaRefs');
    ylabel('Beta');
    xlabel('Image');
    
    subplot(235);
    bar(betaVals,sum(betaRefs,2));
    axis([0 180 0 inf]);
    ylabel('Frequency');
    xlabel('Beta');
    %%
    
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
    subplot(236);
    imags(gammas,betas,avgP.^.2);
    xlabel('Gamma');
    ylabel('Beta');
    colormap jet
%     colorbar;

