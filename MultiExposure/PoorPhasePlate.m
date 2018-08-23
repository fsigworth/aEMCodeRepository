% PoorPhasePlate.m
% Demonstrate phase-plate transfer functions



fmax=0.4;
df=1e-4;
freqs=(0:df:fmax)';

cpars.lambda=EWavelength(200);
cpars.defocus=1;
cpars.Cs=2;
cpars.B=50;
cpars.alpha=.07;
cd=ContrastTransfer(freqs,cpars);


f0=.002;
df=f0/4;
cparsp=cpars;
cparsp.defocus=0;
cparsp.alpha=pi/2;
cparsp.Cs=7;
cp=ContrastTransfer(freqs,cparsp);
cp=cp.*(.5+.5*erf((freqs-f0)/df));

plot(sqrt(freqs),abs([(cd) (cp)]));
SquareRootTicks(freqs);
ylabel('Contrast transfer function');
xlabel('Spatial frequency, Å^{-1} (sqrt scale)');
legend('Defocus contrast, 1 µm','Phase-plate 50nm cut-on');

figure(2);
cdf=cd;
[v x]=min(cdf);
cdf(1:x)=v;
h=cdf./cd;
h(x:numel(h))=1;
sfreqs=sqrt(freqs);
subplot(2,1,1);
plot(sfreqs,h);
ylabel('Filter transfer function');
SquareRootTicks(freqs);
subplot(2,1,2);
plot(sfreqs,[abs(cd) abs(cdf)+.005]);
ylabel('Contrast transfer function');
SquareRootTicks(freqs);
xlabel('Spatial frequency, Å^{-1} (sqrt scale)');
legend('Normal CTF','Boosted');
