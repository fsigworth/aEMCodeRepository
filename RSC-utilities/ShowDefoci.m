% look up defocus values and make histograms

d1=[]; 
d2=[];
nmi=numel(mis);
% for i=1:816; 
for i=1:nmi 
    mi=mis{i}; 
    if isfield(mi,'ctf')
        d1(end+1)=mi.ctf(1).defocus; 
        d2(end+1)=mi.ctf(2).defocus; 
    end; 
end;
bins1=0:.1:4;
figure(5);
subplot(211);
h1=hist(d1,bins1);
subplot(212);
h2=hist(d2,bins1+8);

figure(5);
plot(bins1,[cumsum(h1); cumsum(h2)]')