% k2MapDefects

cd('/Volumes/D215/160909/KvLipo121_2/movie_frames/sq02_1');
load defectLocs.mat
m=ReadEMFile('CountRef_Sep09_17.23.21.dm4');
n=size(m);
%%
def=zeros(n);
locs=round(defectLocs/2);
nanLocs=any(isnan(locs),2);
locs(nanLocs,:)=[];
for i=1:size(locs,1)
def(locs(i,1),locs(i,2))=1;
end;

%%
ndis=64;
m(:)'*def(:)
cc1=Crop(fftshift(real(ifftn(fftn(m).*conj(fftn(def))))),ndis);
imags(cc1);
title('cc1');
drawnow;
% 
% rm=flipud(m);
% rm(:)'*def(:)
% cc2=Crop(fftshift(real(ifftn(fftn(rm).*conj(fftn(def))))),ndis);
% imags(cc2);
% title('cc2');
% 
% drawnow;
% 
% rm=fliplr(rm);
% rm(:)'*def(:)
% cc3=Crop(fftshift(real(ifftn(fftn(rm).*conj(fftn(def))))),ndis);
% imags(cc3);
% title('cc3');
% drawnow;
% 
% rm=fliplr(m);
% rm(:)'*def(:)
% cc4=Crop(fftshift(real(ifftn(fftn(rm).*conj(fftn(def))))),ndis);
% imags(cc4);
% title('cc4');
% drawnow;

mv=ReadMovie('Sep09_17.23.21.tif');
%%
mvs=sum(single(mv),3);
mvs=rot90(mvs,3);

mv2=ReadMovie('Sep09_19.05.10.tif');
mvs2=sum(single(mv2),3);
mvs2=rot90(mvs2,3);
