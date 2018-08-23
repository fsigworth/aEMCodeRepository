% f2MakeFakeMovie
% Create a movie file with known x-shifts.
name='Fakez.tif';
sigmaN=1;
amp=.3;
pixA=1.4;
n=4096;
nim=7;
sx=[7 5 3 2 1 1 0];
noise=sigmaN*single(randn(n,n,nim));
c=CTF(n,pixA,EWavelength(200),2,2,40,.02);
signal=amp*real(ifft2(fft2(randn(n,n))...
    .*ifftshift(c)));
m=zeros(n,n,nim,'single');
for i=1:nim
    m(:,:,i)=circshift(signal,[sx(i),0]);
end;
m=m+noise+1;
disp('Writing');
WriteZTiff(m,pixA,name);
disp('done.');
