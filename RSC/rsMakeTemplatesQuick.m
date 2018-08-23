function templates=rsMakeTemplatesQuick(templateAngles,map)
% This is negligibly faster than rsMakeTemplates, so isn't worth the
% complication.

n=size(map,1);
nangs=size(templateAngles,1);
ks=3;
templates=single(zeros(n,n,nangs));

betaGamma=templateAngles(:,2)+1000*templateAngles(:,3);
[bgSort indSort]=sort(betaGamma);
sortAngles=templateAngles(indSort,:);

bgSteps=find([-inf ; bgSort] ~= [bgSort ; inf]);
nAlphas=diff(bgSteps);
bgSteps(numel(bgSteps))=[];

comp=gridMakePreComp(n,ks);  % Make the pre-compensation function (a 1D array)
F3=gridMakePaddedFT(map,'grid',comp);  % get the 3D fft in a form for slicing.

nSteps=numel(bgSteps);
temps=single(zeros(n,n,nSteps));
for i=1:nSteps
    ind=bgSteps(i);
    angs=rsDegToEuler([0 sortAngles(ind,2:3)]); % skip alpha
    P2=gridExtractPlaneE(F3,angs,ks);  % angs is a 3x1 vector of Euler angles (radians)
    im1=gridRecoverRealImage(P2);
    temps(:,:,i)=im1;
end;


for i=1:nSteps
    ind=bgSteps(i);
    alphas=sortAngles((ind:ind+nAlphas(i)-1)',1);
    im1=temps(:,:,i);
    templates(:,:,indSort((ind:ind+nAlphas(i)-1)'))=...
            rsRotateImage(im1,alphas);
end;
