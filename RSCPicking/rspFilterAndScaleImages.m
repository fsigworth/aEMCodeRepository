function [imgs,mulr,addr]=rspFilterAndScaleImages(mi,dis,rscc)
% Load the first two images of the imgs stack.
% if rscc.mVes isn't present, return only the first image.
doVesImage=isfield(rscc,'mVes') && sum(rscc.mVes(:))~=0;
% ndis=dis.ndis;
effWeights=mi.weights;
if sum(effWeights)==0
    effWeights=1;
end;
if dis.filter(4)
    effWeights(2)=dis.filter(4);
end;
ctComp=dis.filter(3)/100; % percent compensation for CTF
fc=mi.pixA*dis.ds./dis.filter; % Convert from A to A^-1
fc(dis.filter==0)=0;           % 0 A -> 0 A^-1
imgs=single(zeros([dis.ndis 1+doVesImage],'single'));
mu=rspCTFInverseFilter(rscc.m0,mi,ctComp,effWeights);
mf=GaussFilt(GaussHP(mu,fc(1)),fc(2));
midpoint=median(mf(:));
range=[midpoint-(dis.contrast(1)/10) midpoint+dis.contrast(2)/10];
mulr=256/diff(range);
addr=-mulr*range(1);
% [imgs(:,:,1),mulr,addr]=imscale(mf,256,dis.contrast);
imgs(:,:,1)=mf*mulr+addr;
if doVesImage
    mVesf=GaussFilt(rspCTFInverseFilter(rscc.mVes,mi,ctComp,effWeights),fc(2));
    mVesf=GaussHP(mVesf,fc(1));
%     imgs(:,:,2)=imscale(mf-mVesf,256,dis.contrast);
    imgs(:,:,2)=(mf-mVesf)*mulr+addr;
else
    imgs(:,:,2)=0;
end;
