% rsPickerTrends
% Look at how the min amplitude and max spect thresholds vary with defocus.

load('/Users/fred/EMWork/Hideki/161101/KvLipo122_4b/sq10_2w10p192m2tsi.mat')

%%
nmi=numel(si.mi);
nMiParticles=zeros(nmi,1);
miDatOk=false(nmi,1);
pars=zeros(nmi,20);
defs=zeros(nmi,1);
for i=1:nmi
    nParts=sum(si.miIndex==i);
    nMiParticles(i)=nParts;
    mi=si.mi{i};
    particle=mi.particle;
    miDatOk(i)=(nParts>0) && isfield(particle,'autopickPars')...
        && numel(particle.autopickPars)>0;
    if miDatOk(i)
        pars(i,:)=particle.autopickPars;
    end;
    defs(i)=mi.ctf(1).defocus;
end;

amp=pars(miDatOk,1);
amp(amp<1.5)=NaN;
amp(amp>2.5)=NaN;

spect=pars(miDatOk,10);
spect(spect>20)=NaN;
def=defs(miDatOk);

figure(1)
clf;
plot([def amp spect/5]);
legend('defocus','amp','spect/5');

figure(2);
subplot(2,1,1);
plot(def,amp,'o');
xlabel('Defocus');
ylabel('Particle amp');

subplot(2,1,2);
plot(def,spect,'o');
xlabel('Defocus');
ylabel('Spectrum')


% Amp trend: 2.2, 2.1; 3.3,1.65
ampSlope=(1.65-2.1)/(3.3-2.2)
% =-.41
% Spec trend:2.1, 11; 3.5,16
specSlope=(16-11)/(3.5-2.1)
% =3.6
