% reShowRawLProbs.m
% run after reEMReconstruct iteration to look at priors and likelihoods
for iInd=1:1
exf=.1;
exfr=.1;
nt=2*maxShift+1;
nAlphas=ri.angleN(1)*2;
nGammas=ri.angleN(3);
nBetas=ri.angleN(2);

logPX=reshape(rawLProbs.logPXs(:,:,iInd),nt,nt,nAlphas,nGammas,nBetas);
logPrior=reshape(rawLProbs.logPriors(:,:,iInd),nt,nt,nAlphas,nGammas,nBetas);
logPriorG=squeeze(logPrior(:,:,:,1,:));
mxPx=squeeze(max(max(logPX,[],2),[],1));
mxPxA=max(mxPx(:));
mxEPx=exp((mxPx-mxPxA)*exf);

mxPr=squeeze(max(max(logPrior,[],2),[],1));
mxPrA=max(mxPr(:));
mxEPr=exp((mxPr-mxPrA)*exfr);

[val,xiAlpha,xiGamma,xiBeta]=max3d(mxPx);


subplot(3,2,1);
imac(mxEPx(:,xiGamma,:)*256);
xlabel('iAlpha');
ylabel('iBeta');
title(['iGamma = ' num2str(xiGamma)]);

subplot(3,2,2);
imac(mxEPr(:,xiGamma,:)*128+mxEPx(:,xiGamma,:)*128);
title('Prior');

xs=-maxShift:maxShift;
subplot(3,2,3);
xImg=exp((logPX(:,:,xiAlpha,xiGamma,xiBeta)-mxPxA)*exf);
imac(xs,xs,256*xImg);
xlabel('Trans X');
ylabel('Trans Y');
title(['iBeta = ' num2str(xiBeta)]);

subplot(3,2,4);
[val,riAlpha,riGamma,riBeta]=max3d(mxPr);  % is indep. of gamma anyway
riAlpha
riBeta

rImg=exp((logPrior(:,:,xiAlpha,xiGamma,xiBeta)-mxPrA)*exfr);
imac(xs,xs,128*rImg+128*xImg);

logPX(1,1,:)=max(logPX(:));
logPrior(1,1,:)=max(logPrior(:));

subplot(3,2,5);
imacs(imgs(:,:,iInd));

pause(0.2);
end;