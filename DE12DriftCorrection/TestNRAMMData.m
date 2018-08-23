iMean=mean(single(s.imgs),3);
iMean=iMean-mean(iMean(:));

% cc=real(ifftn(fftn(iMean).*conj(fftn(s.br2))));
k=5;
kd=1;
i1=(single(s.imgs(:,:,3))-k*s.dr1)./(s.br1-kd*s.dr1);
% i1=single(s.imgs(:,:,3));
% i1=s.br1-0*s.dr1;
i1(2287,:)=mean(i1(2286:2:2288,:));
i2=(single(s.imgs(:,:,6))-k*s.dr2)./(s.br2-kd*s.dr2);
% i2=single(s.imgs(:,:,6));
% i2=s.br2-0*s.dr2;
i2(2287,:)=mean(i2(2286:2:2288,:));
i1=RemoveOutliers(i1);
i1=i1-mean(i1(:));
i2=RemoveOutliers(i2);
i2=i2-mean(i2(:));

csp=fftn(i1).*conj(fftn(i2));
cc=real(ifftn(csp));
cc1=Crop(fftshift(cc),256);
subplot(2,2,2);
imacs(abs(cc1));
subplot(2,2,4);
plot(sect(cc1));
%%
dr1=RemoveOutliers(s.dr1);
br1=RemoveOutliers(s.br1);
%%
q=(iMean-dr1)./(br1-dr1);
qr=RemoveOutliers(q);
sd=std(qr(:));
md=median(qr(:));
qr(abs(qr-md)>5*sd)=md;
subplot(1,2,1);
imacs(GaussFilt(qr',.05,1)');
subplot(1,2,2);
plot(GaussFiltDCT(median(qr,2),.05));

%% Anisotropic filter
