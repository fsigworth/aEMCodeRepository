% miPickingStats.m
% Plot statistics from SimpleRSPicker robo-picking with generous
% thresholds.

figure(2);
clf;
% Assume we've run MiLoadAll.m
nmi=numel(allMis);
defs=zeros(nmi,1);
nVes=zeros(nmi,1);
for i=1:nmi
    defs(i)=allMis{i}.ctf.defocus;
    nVes(i)=numel(allMis{i}.vesicle.x);
end;

amps=zeros(0,1);
spects=zeros(0,1);
medAmps=zeros(nmi,1);
medSpects=zeros(nmi,1);
for i=5:5
    picks=allMis{i}.particle.picks;
    if size(picks,2)>=8
    amps=[amps; picks(:,5)];
    medAmps(i)=median(picks(:,5));
    spects=[spects; picks(:,8)];
    medSpects(i)=median(picks(:,8));
    end;
end;
ok=amps>0;
amps=amps(ok);
% lamps=log10(amps);
spects=spects(ok);
x0=min(amps);
x1=max(amps);
% lx0=min(lamps);
% lx1=max(lamps);
y1=max(spects);
mysubplot(2,2,1);
plot(amps,spects,'.');
axis([x0 x1 0 y1]);

nBins=ceil(4*sqrt(numel(amps)));
xBins=x0:(x1-x0)/(nBins-1):x1;
yBins=0:y1/(nBins-1):y1;

b=0;
amin=x0;
smax=y1;
while b~='q'
    ampsOk=amps>=amin;
    spectsOk=spects<=smax;
    mysubplot(2,2,2);
    hs=hist(spects(ampsOk),yBins);
    barh(yBins,sqrt(hs));
    mysubplot(2,2,3);
    ha=hist(amps(spectsOk),xBins);
    bar(xBins,sqrt(ha));
    [amin,smax,b]=ginput(1);
end;
totalParticles=sum(ampsOk&spectsOk)