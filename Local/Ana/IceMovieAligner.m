% IceMovieAligner
nc=32;

load('/Users/fred/Box Sync/VesicleMovies/141120View/movie_frames/mb4.mat');
m=single(mb4);
[n,n,ni]=size(m);
mc=m-mean(m(:));
mSum=sum(mc,3);

fm=fft2(mc);
ccs=zeros(nc,nc,ni);
fms=fm;
mca=mc;
ctr=n/2+1;
%%
for i=2:ni
    cc=fftshift(real(ifftn(fm(:,:,i).*conj(fm(:,:,1)))));
    ccs(:,:,i)=Crop(cc,nc);
    [~,ix,iy]=max2di(cc);
    shift=[ix,iy]-ctr;
    fms(:,:,i)=fm(:,:,i).*FourierShift(n,-shift);
    mca(:,:,i)=real(ifftn(fms(:,:,i)));
end;

save('/Users/fred/Box Sync/VesicleMovies/141120View/movie_frames/mca.mat','mca');

mAli=real(ifftn(sum(fms,3)));
figure(1);
imags(mAli)
%% Look at the mean intensity in two spots.
ns=2;
spots=zeros(ns,ni);
spots(1,:)=mean(mean(mca(121:170,735:770,:)));  % broken carbon
spots(2,:)=mean(mean(mca(351:550,491:625,:)));  % central hole
spots(3,:)=mean(mean(mca(351:360,491:500,:)));  % part of central hole

figure(2);
plot(spots');

%%
figure(3);
imovie(mca,.2);

%%  Subtract the mean and look for temporal patterns
mcDiff=mca;
for i=1:ni
    mcDiff(:,:,i)=mca(:,:,i)-mAli/ni;
end;
msk=GaussFilt(mAli,.02)>900;
msk=(msk-GaussFilt(msk,.1))>.1;
for i=1:ni
imaga((1-msk)*128+10*GaussFilt(mcDiff(:,:,i),.05));
pause(0.1);
end;
