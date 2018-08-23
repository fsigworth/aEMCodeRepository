% CC1d
% Look at errors in locating the peak of the cross-correlation function.

n=512;
a=160;
d=40;
sigma=300;
pixA=1.25;
nTrials=10000;
f0=0.01;
c=CTF(n,pixA,.025,4,2,300,.05);
v0=VesicleDensity(n,a,d);


vf=real(ifftn(fftn(v0).*ifftshift(c)));  % ctf-filtered vesicle
vfn=vf+sigma*randn(n);  % A representative noisy image

xs=-n/2:n/2-1;
xctr=n/2+1;

errors=zeros(nTrials,1);
peaks=zeros(nTrials,1);
bins=(-n/2:n/2-1)';
freqs=bins/n;

if f0~=0
    filter=ifftshift(1./(1+(freqs/f0).^2));
else
    filter=1;
end;

vf1=sect(vf);
varvf1=vf1'*vf1;

for i=1:nTrials
    y=vf1+sigma*randn(n,1);
    ref=vf1;
    cc=fftshift(real(ifft(fft(y).*filter.*conj(fft(ref)))))/varvf1;
    [peakVal,peakPos]=max(cc);
    peakPos=peakPos-n/2-1;
    peaks(i)=peakVal;
    errors(i)=peakPos;
end;

subplot(2,2,2);
h=hist(errors,bins);
hmax=max(h);
bar(bins,h);
p1=find(h>hmax/2,1);
p2=find(h>hmax/2,1,'last');
width=p2-p1+1;
title(['Peak: ' num2str(mean(peaks)) '   width at 1/2:  ' num2str(width)]);

subplot(2,2,3);
imags(vfn)

subplot(2,2,4);
plot([y ref]); % Show an example of zero-noise and noisy traces
subplot(2,2,1);
plot(xs,cc,'.-',0,cc(xctr),'r+');
title(peakPos)
drawnow;
