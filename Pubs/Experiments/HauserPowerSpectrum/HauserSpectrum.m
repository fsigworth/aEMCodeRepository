m=imread('~/Desktop/Screen Shot 2019-03-04 at 5.45.31 PM.png');
ms=sum(single(m),3);
ms=rot90(ms,3);
msc=ms(701:1500,101:900);
msc=min(msc,520);

q=fftn(msc);
q(1)=0;
qs=fftshift(q);
% 1 pixel = 10nm
fs=(0:399)/(800*.01); % per um
sq1=Radial(abs(qs).^2);

subplot(221); imags(ms); axis equal
title('Hauser et al. Fig 2C');

subplot(222); imags(msc); axis equal

subplot(223);
imags(Crop(abs(qs).^.4,400));

subplot(224);
semilogy(fs(1:200),sq1(1:200))
hold on;
x=1000/180; % 180 nm periodicity.
plot([x x],[1e8 1e14]);
text(x,1e9,'180nm^{-1}');
hold off;
xlabel('Spatial frequency, um^{-1}');
ylabel('Spectral density');
subplot(221);

