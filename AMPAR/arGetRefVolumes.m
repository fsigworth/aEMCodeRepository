function [vols, mbnOffsetA, volSD]=arGetRefVolumes(pixA,n,mapName,nVols,zSqueeze)
%  Get AMPAR reference volumes from the crystal structure
%  We read the map file stored in the same directory as this function.
%  mbnOffset is the distance of the membrane center from the center of the
%  volume (negative in the case of the AMPAR).
defaultMbnOffsetA=0;

if nargin<4
    nVols=1;
end;
if nargin<5
    zSqueeze=1.2;
end;

mbnHalfThickness=30;

n=n(1);  % volumes are cubes.
pa=AddSlash(fileparts(which('arGetRefVolumes'))); % Get our directory
[pa2,nm,ext]=fileparts(mapName);
if strcmp(ext,'.mat')
    map=load([pa mapName]);  % load the 'map','pixA' and 'mbnOffsetA' variables
    mbnOffsetA=map.mbnOffsetA;
    pixA0=map.pixA;
    map=map.map;
else
    [map,pixA0]=ReadEMFile([pa mapName]);
    mbnOffsetA=defaultMbnOffsetA;
    disp(['Setting the membrane offset to ' num2str(mbnOffsetA)]);
end;

% Set up to squeeze the structure in Z to make more volumes
nctr=ceil((n+1)/2);
nPad=n*1.5;
mag=pixA0/pixA;  % magnification factor
padMap=DownsampleGeneral(map,nPad,mag)/mag; % scale amplitudes by mag also.
mbnOffset=mbnOffsetA/pixA;
mbnTop=round(mbnOffset+mbnHalfThickness/pixA);  % pixel position of top of membrane
mTop=mbnTop+nctr;  % absolute pixel position of membrane top.

v=Crop(padMap,n);
volSD=zeros(nVols,1);
volSD(1)=sqrt(v(:)'*v(:));
vols=zeros([n n n nVols],'single');
vols(:,:,:,1)=v;

xyExpand=sqrt(zSqueeze);
tempVol=v;
tempVol(:,:,n+1:ceil(n*zSqueeze))=0;
% Compress the structure above the membrane by zSqueeze, and expand that
% region horizontally by xyExpand.
for i=2:nVols
  for z=mTop+1:n
      zSamp=mTop+round((z-mTop)*zSqueeze); % We just do nearest-neighbor in z
      v(:,:,z)=DownsampleGeneral(tempVol(:,:,zSamp),n,xyExpand);
  end;
  vols(:,:,:,i)=v;
  tempVol=v;
  tempVol(:,:,n+1:ceil(n*zSqueeze))=0;
  volSD(i)=sqrt(v(:)'*v(:)/numel(v));
end;
