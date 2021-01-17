% rsPickerTrends
% Look at how the min amplitude and max spect thresholds vary with defocus.

% load('/Users/fred/EMWork/Hideki/161101/KvLipo122_4b/sq10_2w10p192m2tsi.mat')
si.mi=allMis;
%%
nmi=numel(si.mi);
nMiParticles=zeros(nmi,1);
miDatOk=false(nmi,1);
pars=zeros(nmi,20);
defs=zeros(nmi,1);
for i=1:nmi
    mi=si.mi{i};
    if ~isfield(mi,'particle') || ~isfield(mi.particle,'picks') || size(mi.particle.picks,1)==0
        continue;
    end;
    flags=mi.particle.picks(:,3);
    okParts= flags==16 | flags==32;
    nParts=sum(okParts);
    nMiParticles(i)=nParts;
    particle=mi.particle;
    miDatOk(i)= (nParts>0) && isfield(particle,'autopickPars')...
        && numel(particle.autopickPars)>0;
    if miDatOk(i)
        pars(i,:)=particle.autopickPars;
    end;
    defs(i)=mi.ctf(1).defocus;
end;


def=defs(miDatOk);
amps=pars(miDatOk,1);
spect=pars(miDatOk,10);
figure(1);
plot([amps spect def]);

ampsOk=amps>.25;
amps=amps(ampsOk);
ampDefs=def(ampsOk);

% amp(amp>1)=NaN;

spectOk=spect<.37;
spects=spect(spectOk);
spectDefs=def(spectOk);

% Do least squares on amps
d0=ampDefs*0+1;
d1=ampDefs-2;
d2=d1.^2;

F=[d0 d1 d2];
ampA=LinLeastSquares(F,amps);
ac=ampA(2:end)/ampA(1)
fitA=F*ampA;
ac=MyInput('Amp coeffs ',ac')';
modFitA=F*(ampA(1)*[1 ac']');

figure(2);
subplot(2,1,1);
plot(ampDefs,amps,'o');
hold on; plot(ampDefs,[fitA modFitA],'.');
hold off;
xlabel('Defocus');
ylabel('Particle amp');


d0=spectDefs*0+1;
d1=(spectDefs-2);
d2=d1.^2;

F=[d0 d1 d2];
spectA=LinLeastSquares(F,spects);
sc=spectA(2:end)/spectA(1)
fitS=F*spectA;
sc=MyInput('Spect coeffs ',sc')';
modFitS=F*(fitS(1)*[1 sc']');

subplot(2,1,2);
plot(spectDefs,spects,'o');
hold on
plot(spectDefs,[fitS modFitS],'.');
hold off;
xlabel('Defocus');
ylabel('Spectrum')


% Amp trend: 2.2, 2.1; 3.3,1.65
ampSlope=(.4 - .33)/(1.5-3.30)
% ampSlope=(1.65-2.1)/(3.3-2.2)
% =-.41
% Spec trend:2.1, 11; 3.5,16
% Spec trand .18,1.5; .26,3
specSlope=(.18-.26)/(1.5-3)
% =3.6
