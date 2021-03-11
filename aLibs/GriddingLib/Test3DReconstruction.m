% Test3DReconstruction.m
% Test the gridding-extraction and insertion functions
% by doing a random tomography simulation.
%

n=128;   % number of points
na=100;  % number of angles

% Make a phantom
w=2;  % sharpness of objects
origVol=fuzzymask(n,3,n*.25,w,n*[.4 .4 .5])+2*fuzzymask(n,3,n*.1,w,n*[.5 .8 .5]);

figure(1);

ShowSections(origVol);
drawnow;

% Make projections
comp=gridMakePreComp(n,3);
P0=gridMakePaddedFT(origVol,'grid',comp);
% random angles
angles=[2*pi*rand(na,1) acos(rand(na,1)) 2*pi*rand(na,1)];

disp('Making projections');
tic
imgs=zeros(n,n,na,'single');
for i=1:na
    p=gridExtractPlane(P0,angles(i,:));
    imgs(:,:,i)=gridRecoverRealImage(p);
end;
toc

figure(2);
ImagicDisplay2(BinImage(imgs,4));
drawnow;

disp('Reconstruction');
tic
pc=gridMakeConstFT(n,2,1);  % const 2D plane
P1=gridMakeConstFT(n,3,0);  % 3D Fourier volume
PN=P1;                      % 3D normalization volume

for i=1:na
    p=gridMakePaddedFT(imgs(:,:,i));
    P1=gridInsertPlane(p,P1,angles(i,:),3);
    PN=gridInsertPlane(pc,PN,angles(i,:),3);
end;
%%
figure(3);  % Show the Fourier space coverage
ShowSections(PN.PadFT.*Radius3(PN.np1));

%%
k=.001;  % wiener constant for Fourier normalization
PR=P1;
PR.PadFT=P1.PadFT.*PN.PadFT./(k+PN.PadFT.^2);
comp=gridMakePreComp(n,3);
vol=gridRecoverRealImage(PR,comp);
toc

figure(4);
clf;
ShowSections(vol);
