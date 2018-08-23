imgSize=64;
model=zeros([imgSize imgSize]);
[x,y]=meshgrid(1:imgSize,1:imgSize);
r=sqrt((x-imgSize/2-0.5).^2+(y-imgSize/2-0.5).^2);
model=double((r>20)&(r<24));
f1=fspecial('gaussian',11,0.8);
f2=fspecial('gaussian',11,2);
f=f1-0.95*f2;
model=conv2(model,f,'same');
%imshow(mat2gray(model))
rings=makeRings(imgSize,12);

sigmas=1;
for iSigma=1:numel(sigmas);
%Create the two images
sigma=sigmas(iSigma);
 img1=model+sigma*randn(size(model));
 figure(10);
 imshow(mat2gray(img1));
 figure(2);
 
% img2=model+sigma*randn(size(model));
% img2=circshift(img2,10,1);
sIndex=[];
fIndex=[];
bIndex=[];
maxShift=16;
for  i=1:1000
    img1=model+sigma*randn(size(model));
    img2=model+sigma*randn(size(model));
    img2=circshift(img2,2,1);
    fIndex=[fIndex frcMatch(img1,img2,maxShift)];
    sIndex=[sIndex shiftMatch(img1,img2,maxShift)];
    bIndex=[bIndex bayesShiftMatch(model,img2,maxShift,2*sigma)];
end

% figure(2);
% plot(fIndex,sIndex,'*');
%figure(3);
hf=hist(fIndex,[1:maxShift]);
%figure(4)
hs=hist(sIndex,[1:maxShift]);
hb=hist(bIndex,[1:maxShift]);
figure(3);
subplot(2,2,iSigma)
hold off;
plot([1:maxShift],hf,'b*-');
hold on;
plot([1:maxShift],hs,'r*-');
plot([1:maxShift],hb,'g*-');

legend({'frc' 'shift' 'bayes'});
title(['\sigma = ' num2str(sigma)]);
end;

