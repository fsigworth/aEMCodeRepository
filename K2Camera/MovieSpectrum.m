% MovieSpectrum

flagSegments=[1 20; 23 inf];
    [fname, pa]=uigetfile('*.mrc','Select mrc file');
    if isnumeric(pa) % File selection cancelled
        return
    end;
    cd(pa);
                m=single(ReadMRC(fname));
disp('Removing outlers');
thresh=Percentile(fst,1-1e-7);
fst(fst>thresh)=thresh;
[nx,ny,nim]=size(m);
n0=NextNiceNumber(max(nx,ny));
% Make the padded stack mc
mc=zeros(n0,n0,nim,'single');
for i=1:nim
    q=m(:,:,i);
    me=mean(q(:));
    mc(:,:,i)=Crop(q-me,n0);  % subtract mean
end;

fmc=fft2(mc);  % Fouriers
nw=n0/ds;
fmd=Crop(mc,nw,1);  % crop the whole stack

fmask=fuzzymask(nw,2,ro,ro/10);


%%
disp('Binning');
for i=1:nim
    q=single(st(:,:,i));
    me=mean(q(:));
    stbin(:,:,i)=SharpFilt(BinImage(Crop(q,n0,0,me),ds1),.5/ds2,.1/ds2);
end;



m=Crop(m,n0,1);
m=min(m,6);


mb=BinImage(m,nbin);  % bin the stack




%%
% Real-space correlations
mb=mb-mean(mb(:));
[nbx,nby,nim]=size(mb);
mbsum1=zeros([nbx nby nim]);
mbsum2=mbsum1;
sum1=zeros([nbx nby]);
sum2=zeros([nbx nby]);

for i=1:nim
    j=num-i+1;
    mbsum1(:,:,i)=sum1;
    mbsum2(:,:,j)=sum2;
    fm1=GaussHP(mb(:,:,i),fch);
    fm2=GaussHP(mb(:,:,j),fch);
    corr1=fm1(:)'*sum1(:);
    corr2=fm2(:)'*sum2(:);
    
    
    sum1=sum1+mb(:,:,i);
    sum2=sum2+mn(:,:,j);

    
%%
j=2;
m=m-mean(m(:));
p1=flagSegments(j,1);
p2=min(nim,flagSegments(j,2));

fm=fftshift(ffts(m));




sp1=fftshift(abs(fft2(sum(m(:,:,p1:p2),3)).^2));
sp2=fftshift(sum(abs(fft2(m(:,:,p1:p2))).^2,3));





%%
indivSpectra=abs(fft2(m)).^2;
%%
littleSpectra=zeros(floor([nx/nbin ny/nbin nim]));
sumSpectra=zeros(floor([nx ny]/nbin));
msk=1-fuzzymask([nx ny],2,20,1);
msk2=1-fuzzymask([nx ny],2,nx/4,1);
mmean=mean(m(:));
for i=1:nim
    sp=BinImage(msk.*indivSpectra(:,:,i)+mmean*(1-msk),nbin);
    littleSpectra(:,:,i)=sp;
    sumSpectra=sumSpectra+sp;
    if i>1
    corr(i)=sumSpectra(:)'*sp(:)/(i-1);
    end;
end;
plot(corr);

return



figure(1);
SetGrayscale;
ndis=[nx ny]/2;
nbin=4;
dexp=.2;

msk=1-fuzzymask([nx ny],2,20,1);
q=sp1.*msk;
spmean=sum(q(:))/sum(msk(:));
sp1=sp1.*msk+spmean*(1-msk);
subplot(221);

sp1b=BinImage(Crop(sp1,ndis),nbin);
imacs(sp1b.^dexp);
subplot(223);
plot(Radial(sp1b));

q=sp2.*msk;
spmean=sum(q(:))/sum(msk(:));
sp2=sp2.*msk+spmean*(1-msk);
subplot(222);
sp2b=BinImage(Crop(sp2,ndis),nbin);
imacs(sp2b.^dexp);
subplot(224);
plot(Radial(sp2b));