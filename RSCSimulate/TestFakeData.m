% TestFakeData

load /Volumes/TetraData/Structures/AMPAR/3KG2RotMap5.8A.mat
nt=size(map,1);
map=map/64;  % approx scaling.

membraneOffset=-24/2;  % downsampled map by 2.
dGamma=22.5;
n=256;
r0=80;
vesAmp=.01;
% radius is r0*5.8A.
ctr=(n/2+1)*[1 1];

particleAngs=[0 70 0];
particleAngs=[180 90 30; 90 40 00; -90 90 0];
particleAngs=[90 40 00; 0 0 0; -90 90 0];
rso=[1 1 1];
orientAngs=particleAngs;
for i=1:size(particleAngs,1)
    if rso(i)<0  % inside out
        angs=particleAngs(i,:);
        orientAngs(i,:)=[180+angs(1) 180-angs(2) angs(3)];
    end;
end;
    
templates=rsMakeTemplates(orientAngs,map);
img=zeros(n,n);
img=vesAmp*VesicleDensity(n,r0,10);
nparts=size(particleAngs,1);
shifts=zeros(nparts,2);
iGammas=round(orientAngs(:,3)/dGamma)+1;
iGammas(iGammas>8)=0;
for i=1:size(particleAngs,1)
    angs=particleAngs(i,:);
    shift=(r0-membraneOffset*rso(i))*sind(angs(2))*[sind(angs(1)) cosd(angs(1))]...
        +ctr;
shifts(i,:)=shift;
particleImgs=ifftshift(Crop(templates(:,:,i),n));
    img=img+real(ifftn(fftn(particleImgs).*FourierShift(n,shift)));
end;
figure(2);
img(37,:)=1;
imacs(img);
