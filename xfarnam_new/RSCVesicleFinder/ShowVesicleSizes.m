% ShowVesicleSizes

% if nargin<1  % put up a file selector
[fname pa]=uigetfile('*mi.mat','Select an mi file');
% pa
if ~ischar(fname)
    return
end;

    load([pa '/' fname]);

%%
minS=.03;
maxS=.042;
maxR=500;
nbins=30;

diams=2*mi.vesicle.r*mi.pixA;
subplot(2,1,1);
set(gca,'fontsize',14);
plot(diams, mi.vesicle.s,'.','markersize',10);
axis([-inf maxR minS maxS]);
title(fname,'interpreter','none');
xlabel('Vesicle diameter, Å ');
ylabel('Amplitude ');


subplot(2,1,2);
set(gca,'fontsize',14);
selDiams=diams(mi.vesicle.s>minS & mi.vesicle.s<maxS);
hist(selDiams,nbins);
xlabel('Vesicle diameter, Å ');
ylabel('Frequency');
