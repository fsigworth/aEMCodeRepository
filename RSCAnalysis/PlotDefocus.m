% PlotDefocus.m

load Info/allMis.mat
nmi=numel(allMis);
d=zeros(nmi,1);
for i=1:nmi
    mi=allMis{i};
%     if isfield(mi,'ctf') && isfield(mi.ctf,'defocus')
%         d(i)=mi.ctf(1).defocus;
%     end;
d(i)=mi.doses(1);
end;
%%
subplot(211);
plot(d);
ylabel('Defocus, um');
ylabel('Dose');
title(pwd);
subplot(212);
hist(d,50);
ylabel('Defocus, um');
ylabel('Dose');
title([num2str(nmi) ' images']);
