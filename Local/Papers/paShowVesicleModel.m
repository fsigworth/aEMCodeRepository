% paShowVesicleModel.m
pad=15;

vesModel=mi.vesicleModel;
nm=numel(vesModel);
nm=nm+2*pad;
vesModel=[zeros(pad,1,'single'); vesModel; zeros(pad,1,'single')];
pixA=mi.pixA;
nm2=(nm-1)/2;
pVals=(-nm2:nm2)';
zVals=pVals*pixA;

nPks=numel(mi.vesicle.extraPeaks);
pFcns=zeros(nm,nPks);
for j=1:nPks
    pFcns(:,j)=exp(-(pVals-mi.vesicle.extraPeaks(j)).^2/(2*mi.vesicle.extraSD)^2);
end;

plot(zVals,vesModel,'k-','linewidth',1.5);
hold on;
plot(zVals,1.9*pFcns, 'linewidth', 1);
hold off;
axis([-60 60 0 2.2]);
xlabel('Membrane radial distance, Å');
ylabel('Inner potential, V');
