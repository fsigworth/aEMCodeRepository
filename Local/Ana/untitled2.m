% VesicleTracking.m

cd('/Users/fred/Box Sync/VesicleMovies')

load modelVesicle.mat
load Vesicle23.mat

n=256;
ves=Crop(Downsample(modelVesicle,384),n)*8;
%   scale up by 8 allows good subtraction
mvd=Downsample(mvx,384,1);
mvd=Crop(mvd,n,1);
nf=17;  % number of frames
mvd=mvd(:,:,1:nf);  % just the low-defocus part.
mx=sum(mvd,3);

cc=zeros(n,n,nf);
for i=1:nf
    cc(:,:,i)=fftshift(ifftn(fftn(mvd(:,:,i)).*conj(fftn(ves))));
    [mxv,xs(i),ys(i)]=max2di(cc(:,:,i));
    