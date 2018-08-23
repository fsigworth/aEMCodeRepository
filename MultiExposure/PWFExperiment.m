% PWFExperiment

doHuffman=1;

imagePath='/Volumes/TetraData/EMWork/Scott/';
imageName='10sep19a_a_00006gr_00022sq_v01_00002hl_v01_00003en.mrc';
imageName='11oct20b_20170729_04_4096x4096_bright_0.mrc';
floodName='floodbeam.mrc';

name=[imagePath imageName];
% name=[imagePath floodName];
m=double(ReadMRC(name));
me=mean(m(:));
% m=m-me);
% m=mePreWhiten(m);
n=size(m,1);

%%
f=double(Radius(n))/n;
sp2=CCDModelSpectrum(f,1);
sp2=sp2/(sp2(n/2+1,n/2+1));  % normalize to 1
mw=real(ifftn(fftn(m)./ifftshift(sqrt(sp2+.0))));  % pre-whitened
% %
% %
% % sp1=Radial(sp2);
% %
% sp0=fftshift(abs(fftn(mw)).^2);
% imacs(sp0);
% sp1=Radial(sp0);
% % RadialPowerSpectrum(mw);
% figure(1);
% subplot(2,1,1);
% semilogy([sp1]);
% subplot(2,1,2);
% %%
% imacs(abs(GaussFilt(sp0.*(1-fuzzymask(n,2,100,10)),.03)).^.5);
% return

%%
% Try compression
% sensitivity is 33 counts/e; assume quantal sens is 50 counts/e.
infac=8;  % variance factor  % gives 3 pits/pixel, >100:1
infac=4;  % 2 b/pixel, ~25:1
rsens=33*1.5;  % convert to counts
k=2*infac/sqrt(rsens);
msq=sign(mw).* floor(k*sqrt(abs(mw)));
figure(1);
subplot(2,1,1);
hist(msq(:),1000);


% Do the reconstruction
% Find the center of the interval between the squares of
% msk/k (msk/k+1/k)
mss=msq/k;
mrq=sign(mss).*(abs(mss).^2+abs(mss)/k+0.5/(k^2));  % shift to center of bin

subplot(2,1,2);
plot([sect(mw) sect(mrq)]);

plot([sect(mw)-sect(mrq) sect(mw)]);

mr=real(ifftn(fftn(mrq).*ifftshift(sqrt(sp2))));

plot([sect(mr)-sect(m) sect(m)]);


figure(2);
spm=RadialPowerSpectrum(m);
spmr=RadialPowerSpectrum(mr);

semilogy([spm abs(spm-spmr)]);
drawnow;

if doHuffman
    nb=256;  % block size
    blocks=ceil(n/nb);
    
    bR=zeros(blocks,blocks);
    disp('Huffman encoding');
    code=[];
    tic
    for j=1:blocks
        y0=(j-1)*nb+1;
        y1=min(j*nb,n);
        for i=1:blocks
            x0=(i-1)*nb+1;
            x1=min(i*nb,n);
            blockData=msq(x0:x1,y0:y1);
            [codedBlock bitRate]=huff03(blockData(:),0);
            bR(i,j)=bitRate(1);
            disp([i j bitRate(1)]);
            code=[code;uint8(codedBlock)];
        end;
    end;
    toc
    return
    
end;

return
%%

n=4096;

mx=zeros(n,n);
my=mx;

s1=2.5;
s2=-.7;
ct=n/2+1;

mx(ct,ct)=s1+s2+1;
mx(ct+1,ct)=-s1;
mx(ct+2,ct)=-s2;
my(ct,ct)=s+1;
my(ct,ct+1)=-s;

fx=fftn(fftshift(mx));
% fy=fftn(fftshift(my));
sx=sectr(fftshift(fx));
% sz=sectr(fftshift(fx.*fy));
% semilogy(abs([sx sz]));
subplot(2,1,2);
semilogy(abs(sx).^2.*sp0);

