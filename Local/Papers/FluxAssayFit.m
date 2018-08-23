% paFluxAssayFit.m

cd('/Users/fred/aEMCodeRepository/Papers')
load FluxData161129.mat

times=[100 300 500 700]+20;
texts={['<- CCCP  ';'  1uM    '] ['<- Valino  ';'   20 nM   '] ['<- Monensin  ';'      20nM   '] ['<- Monensin  ';'       1uM   ']}

T=vesT;
Y=vesY/1e6;
tText='Empty liposome control';
xMax=800;
yMax=1.8;
yMin=.4;
legLabelY=.8;

% T=kvT;
% Y=kvY(:,1:4)/1e6;
% tText='Kv liposomes';
% xMax=800;
% yMax=6.4;
% yMin=3.7;
% legLabelY=4.44
% 
figure(1);
plot(T,Y,'.-','linewidth',2,'markersize',10);
hold on;
for i=1:numel(times)
    text(times(i),yMax,texts{i},'rotation',90,'horizontalalignment','right');
end;
hold off;
txt=num2str(vesV');
txt=[txt repmat(' mV',4,1)];
legend(txt,'location','southeast');
text(700,legLabelY,'E_K');
ylabel('Flourescence, 10^6 counts');
xlabel('Time, s');
title(tText);
axis([0 xMax yMin yMax]);

kYMax=6e6;
kYmin=4e6;

