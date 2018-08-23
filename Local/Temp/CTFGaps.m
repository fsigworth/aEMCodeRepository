% CTFGaps.m
% Examine the gaps left by CTFs covering a range of defocus values

kV=300;
minD=4;
maxD=5;
nD=20;
spread=.1;
deltaD=(maxD-minD)/(nD-1);
defoci1=(minD:deltaD:maxD)';
nd1=numel(defoci1);
nrep=100;
nd=nrep*nd1;
defoci=repmat(defoci1,nrep,1)+spread*rand(nd,1);
%
% defoci=val;
% nd=numel(val);
%
freqs=(0:.0001:.2)';

lambda=EWavelength(kV);
nf=numel(freqs);

Cs=2;
alpha=.02;
B=0;

cts=zeros(nf,nd);


for i=1:nd
    cts(:,i)=abs(ContrastTransfer(freqs,lambda,defoci(i),Cs,B,alpha));
end;

figure(2);
mysubplot(2,1,1);
plot(freqs,cts);%,'.','markersize',1);
mysubplot(2,1,2);
plot(freqs,mean(cts.^2,2));
xlabel('Spatial frequency, Å^{-1}');
ylabel('Mean squared CTF');
legend([num2str(minD) ' to ' num2str(maxD) 'um;  rand=' num2str(spread)]);