% ComputePhasePlateCTF

cd('/Volumes/D223/170708/KvLipo155_1pp')
allNames=f2FindInfoFiles;

miName=allNames{2};

disp(miName);
mi=ReadMiFile(miName);
n=mi.imageSize;
ds=4;
mScale=1;

%%
m0=meReadMergedImage(mi);
m=Downsample(m0,n/ds)*mScale;
msk=meGetMask(mi,n/ds);
subplot(232);
imags(GaussFilt(m.*msk,.1));
title(miName,'interpreter','none');
drawnow;

% mi.vesicle.s(:,:,2:end)=[];
% mi.vesicle.s(:,2:end)=0;

subplot(236);

disp('Making model vesicles');
v0=meMakeModelVesicles(mi,960,0,0,0);
v=min(v0,0); % force it to be non-positive

imags(v.*msk);
title('Vesicle model');
drawnow;

fm=fftshift(fftn(m.*msk));
fv=fftshift(fftn(v.*msk));

freqs=(0:479)/(mi.pixA*mi.imageSize(1));
%%
figure(1);
disp('Computing 1d spectra');
spFactor=mi.doses(1)/prod(n/ds)*ds^2;
spm=Radial(abs(fm).^2)*spFactor;
spv=Radial(abs(fv).^2)*spFactor;
disp('...done.');
%
subplot(231);
semilogy(freqs,[spm spv]);
xlabel('Frequency, Å^2');
ylabel('Power spectrum');
legend('Image', 'Model');

k=.01;
C=real(fm.*conj(fv))./(k+abs(fv).^2);
% C=real(fvhp.*conj(fv))./(k+abs(fv).^2);
Cf=GaussFilt(C,.1);
c1=Radial(Cf);

subplot(234);
plot(freqs,c1);
xlabel('Frequency, Å^2');
ylabel('Coherence');
%%
subplot(235);
ndis=120;
p=find(c1>.5,1)-1;
fp=freqs(p);
spScale=.1/spm(p);
logSps=min(2,log10([spm(1:ndis) spv(1:ndis)])-log10(spm(ndis)));
c1Fit=.5+.5*erf((freqs(1:ndis)'-fp)/(.8*fp));
plot(freqs(1:ndis),[logSps c1(1:ndis) ],'.-');
hold on;
plot(freqs(1:ndis),c1Fit,'color',[.6 .6 .6]);
plot(freqs(p),c1(p),'k+');
hold off;
legend('Log data spectrum','Log model spectrum','Coherence','Cut-on');
title(['Cut-on ' num2str(freqs(p))]);
xlabel('Frequency, Å^2');
ylabel('Coherence');

fc=freqs(p)*1.5*ds*mi.pixA;
vhp=GaussHP(v.*msk,fc);
subplot(233);
imags(vhp);
title('Filtered vesicle model');
% fvhp=fftshift(fftn(vhp));
