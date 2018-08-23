% RotAutocorr.m

load /Users/fred/EMWork/Simulations/Kv1.2/TotalDens.mat % 1A / pixel
n=240; % pad to a nice number
msk=fuzzymask(n,3,n*0.45,n*0.1);

rvol=Crop(TotalDens,n);
norm=(rvol(:)'*msk(:))/(msk(:)'*msk(:));
rvol=(rvol-norm).*msk;
figure(1);
ShowSections2(rvol);


fvol=gridMakePaddedFT(rvol); 

ks=5;
comp=gridMakePreComp(n,ks);
gammas=(0:89)';
nim=size(gammas,1);
angles=zeros(nim,3);
angles(:,1)=gammas;
angles(:,2)=pi/2;
slices=zeros(n,n,nim,'single');
for i=1:nim
    p2=gridExtractPlane(fvol,angles(i,:)*pi/180,ks);
    slices(:,:,i)=gridRecoverRealImage(p2,comp);
    imacs(slices(:,:,i));
    drawnow;
end;

%%
% Now compute the acf along rotation
figure(4);
SetGrayscale;

subplot(221);
imacs(sum(rvol,3));
subplot(222);
imacs(sum(rvol,2));

slc=(reshape(slices,n^2,nim))';  % rotation is 1st dimension

acf2=real(ifft(abs(fft(slc)).^2));
acf=fftshift(sum(acf2,2));

subplot(223);
plot(acf,'.-','markersize',10);
title('Rotational ACF');
xlabel('Degrees');

sp2=abs(fft(slc)).^2;
sp=mean(sp2,2)/nim^2;
fs=(0:nim/2-1)/nim;  % frequencies in A^-1
subplot(224);
semilogy(fs,sp(1:nim/2).*(0:nim/2-1)');
title('Rotational Spectrum');
