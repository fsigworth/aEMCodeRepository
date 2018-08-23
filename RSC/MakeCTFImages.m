% MakeCTFImages
ns=48;
pixA=5.8;
lambdaE=EWavelength(200);
SigmaN=1;
load 3KG2RotMap2.9A.mat
map2=Downsample(map,40);

map=Crop(map2/2000,48);  % approx normalization to give s/n=10, pad to 48 pixels.
ri.angleStep=10;  % degrees.  5 degrees would be more realistic
ri.symmetry=2;
ri.nGamma=round(90/ri.angleStep);
angleList=rsMakeTemplateAngles(ri);
templates=rsMakeTemplates(angleList,map);
%c=abs(CTF(ns,pixA,lambdaE,def,2,B,.07));

alphaSteps=4;
alphas=(-alphaSteps:alphaSteps)*ri.angleStep;

trans=[0,0];
img60=rsRotateImage(circshift(rsMakeTemplates([0 60 0],map),trans),alphas);
trans=[2,0];
img91=rsRotateImage(circshift(rsMakeTemplates([0 90 0],map),trans),alphas);
trans=[4,0];
img32=rsRotateImage(circshift(rsMakeTemplates([0 30 0],map),trans),alphas);
nimg=45;
imgs=1*reshape([img32,img60,img91,img91,img91],48,48,nimg);
nimg=45*4;
imgs=2*reshape([imgs templates(:,:,1:45) templates(:,:,46:90) templates(:,:,91:135)],48,48,nimg);

rotImgs=imgs;
c=imgs;
%norm=c;
def=zeros(nimg,1);
for i=1:nimg
    def(i)=rand*2+1;  % defocus between 1 and 3
    B=50+50*def(i);
    c(:,:,i)=abs(CTF(ns,pixA,lambdaE,def(i),2,B,.07)); % CTF
    %norms(:,:,i)=fftshift(real(ifftn(ifftshift(c(:,:,i).^2))));
    %rotImgs(:,:,i)=real(ifftn(fftn(imgs(:,:,i)).*ifftshift(c(:,:,i))));
    rotImgs(:,:,i)=fuzzymask(ns,2,ns*0.45,ns*.05)*SigmaN.*randn(ns,ns)+real(ifftn(fftn(imgs(:,:,i)).*ifftshift(c(:,:,i))));
end;