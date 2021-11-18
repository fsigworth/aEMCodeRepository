% HRPickingFigs.m
% Make the simulated vesicle and micelle micrographs

load HRPicking/Figs/TM_Micelle_maps.mat % mMap0 tMapsh micDens1
nC=176;
tMapsh=Crop(tMapsh,nC);
mMap0=Crop(mMap0,nC);
[compMap,sC]=ReadMRC('HRPicking/compMap.mrc');
cMap=circshift(DownsampleGeneral(compMap,176,sC.pixA),[0 0 21]);

%%
nm=size(tMapsh,1);
ctrm=ceil((nm+1)/2);

model=max(0,micDens1(2:end)-.00191);


n=1024;
ctr=ceil((n+1)/2);
a=200;

ves=VesicleFromModel(n,a,model);
% imags(ves)

% figure(5);
% imags(proj)

angs=[0 90 -90; 0 120 220; 0 150 -45];
projs=rlMakeTemplates(angs,tMapsh);

angs2=[0 120 135];
projs(:,:,end+1)=rlMakeTemplates(angs2,cMap);
angsR=[angs;angs2];
img=zeros(n);
angsR(:,3)=270-angsR(:,3);
% angs=[0 160 270-angs(1,3); 0 90 270-angs(2,3)];
for i=1:size(angsR,1)
    Y=a*sind(angsR(i,2));
    pos=[ctr;ctr]+RotMatrix2(angsR(i,3)*pi/180)*[0;Y];
    img=img+ExtractImage(projs(:,:,i),round(pos),n,1);
end;
%%
% Micelle + TM map
angs3=zeros(4,3);
angs3(:,2)=[90 110 135 180]';
angs3(:,3)=-90;
projs3=rlMakeTemplates(angs3,tMapsh+mMap0);
projs3m=rlMakeTemplates(angs3,mMap0);
ys=[150 50 -80 -240]+50;
img3=zeros(n,n,'single');
mic3=zeros(n,n,'single');

for i=1:size(angs3,1)
    img3=img3+ExtractImage(projs3(:,:,i),[ctr ctr+ys(i)],n,1);
    mic3=mic3+ExtractImage(projs3m(:,:,i),[ctr ctr+ys(i)],n,1);
end;
%
%
%
ctPars.defocus=1.4;
ctPars.lambda=EWavelength(300);
ctPars.alpha=.07;
ctPars.Cs=2.7;
ctPars.B=100;
c=CTF(n,1,ctPars);

cImg=real(ifftn(fftn(img).*ifftshift(c)));
cVes=real(ifftn(fftn(ves).*ifftshift(c)));
cImg3=real(ifftn(fftn(img3).*ifftshift(c)));
cMic3=real(ifftn(fftn(mic3).*ifftshift(c)));

sigma=.6;
fc=.15;
dn=640; % no. of pixels to display
noise=sigma*randn(n);

figure(1);

mysubplot(221);
imags(Crop(img+ves,dn));
axis equal off

mysubplot(222);
imags(Crop(GaussFilt(cImg+noise,fc),dn));
axis equal off

mysubplot(223);
imags(Crop(img-ves,dn));
axis equal off

mysubplot(224);
imags(Crop(GaussFilt(cImg-cVes+noise,fc),dn));
axis equal off;

%
% Now do the same for micelles
%
figure(2)
mysubplot(221);
imags(Crop(img3,dn));
axis equal off

mysubplot(222);
imags(Crop(GaussFilt(cImg3+noise,fc),dn))
axis equal off

mysubplot(223)
imags(Crop(img3-mic3,dn));
axis equal off;

mysubplot(224)
imags(Crop(GaussFilt(cImg3-cMic3+noise,fc),dn))
axis equal off


