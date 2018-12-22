% PlotMicrographPeakValues
f0=.0025;  % lowpass Gauss

load Info/allMis.mat
nmi=numel(allMis);
d=zeros(nmi,2);
for i=1162:nmi
    mi=allMis{i};
imName=mi.imageFilenames{1};
[pa,nm,ex]=fileparts(imName);
isName=['Micrograph_ds/' nm 's' ex];
[ms,s]=ReadMRC(isName);
mf=GaussFilt(ms,s.pixA*f0); % approx 400 A filter.
imags(ms);
title(isName);
drawnow;
pause
d(i,1)=max(mf(:));
d(i,2)=median(mf(:));
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
