% NISMaskMaker

inMapName='MapsForSubtraction/210705/cryosparc_P8_J67_010_volume_map_sharp(1).mrc'; % perrhenate
outMaskName='Masks/Mask1.mrc';
outMapName='Masks/MaskedMap1.mrc';
doWrite=1;


[map,s]=ReadMRC(inMapName);

n=size(map,1) % box size
ctr=floor((n+1)/2); % coordinate center
nDis=128; % little box for display

blobCtr=[10 20 0]; % center of one blob; the other is reflected in X,Y
blobRadius=[15 10 10]; % radius in X, Y, Z
fuzz=5; % soft mask width

% make 2 blobs
blobs2=fuzzymask(n,3,blobRadius,fuzz,ctr+blobCtr) ...
    + fuzzymask(n,3,blobRadius,fuzz,ctr-blobCtr);

blobs2=max(0,min(1,blobs2)); % force values to be strictly 0..1 for relion

figure(1);
ShowSections(Crop(map,nDis)); % Show the original map

figure(2); % Show the mask
ShowSections(Crop(blobs2,nDis));

figure(3); % Show mask * map
ShowSections(Crop(blobs2.*map,nDis));

if doWrite
    WriteMRC(blobs2,s.pixA,outMaskName);
    disp([AddSlash(pwd) outMaskName ' written.']);
    WriteMRC(blobs2.*map,s.pixA,outMapName);
    disp([AddSlash(pwd) outMapName ' written.']);
else
    disp('Not written.');
end;

