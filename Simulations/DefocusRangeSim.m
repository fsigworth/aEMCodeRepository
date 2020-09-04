% DefocusRangeSim.m

pixA=4;
n=1024;
freqs=(0:n-1)/(n*pixA);
kV=300;
lambda=EWavelength(kV);

Cs=2.7;
B=0;
alpha=.05;

ndef=100;
std=.4;
me=1.9;
defoci=randn(100,1)*std+me;




ctfs=zeros(n,ndef);
for i=1:ndef
    ctfs(:,i)=ContrastTransfer(freqs,lambda,defoci(i),Cs,B,alpha);
end;

subplot(311);
plot(freqs,ctfs.^2);

subplot(312);
plot(freqs,mean(ctfs.^2,2));

subplot(313);
hist(defoci,20);
