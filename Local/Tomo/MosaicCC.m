function [shift, overlap]=MosaicCC(fPadImgs,fPadMsks,i,j)

[npx,npy,nz]=size(padImgs);
ctr=[npx npy]/2+1;
cc0=fftshift(real(ifftn(fPadImgs(:,:,i).*conj(fPadImgs(:,:,j)))));
ccm=fftshift(real(ifftn(fPadMsks(:,:,i).*conj(fPadMsks(:,:,j)))));
ccq=cc0./(2e4+ccm);
cc=GaussHP(ccq,.1);
subplot(2,2,2);
imacs(cc);
[val, ix, iy]=max2d(cc);
overlap=ccm(ix,iy);
shift=[ix iy]-ctr;
% im2=circshift(padImgs(:,:,j),shift);
%         subplot(2,2,3);
%         im1=padImgs(:,:,i);
%         imacs(im1);
%         subplot(2,2,4);
%         imacs(im2);
%         subplot(2,2,1);
%         imacs(im2+im1);
%         title([i j]);
%         drawnow;
