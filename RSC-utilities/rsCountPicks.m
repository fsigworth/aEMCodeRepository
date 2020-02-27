% rsCountPicks
% Based on rsCountParticles
% Count all the picks in a directory full of info files. Compare these with
% particles in a star file from a selection operation.

useAllMis=0;
infoDir='Info/';
maxEntries=inf;
stride=1000;
allMisName='Info/allMis1.mat';
if useAllMis
    disp(['Loading ' allMisName]);
    load(allMisName);
    nEntries=numel(allMis);
    nEntries=numel(allMis)-1; %%%%%%%%
else
    d=dir(infoDir);
    nEntries=min(numel(d),maxEntries);
end;

siName='Stack1/Dec20_02.10.04p128ztsi.mat';
load(siName);
starName='Stack1/Select/job008/particles.star';
starName='Stack1/Select/job013/particles.star';

[dNames,data]=ReadStarFile(starName);
q=data{1};
nl=numel(q.rlnImageName);
for i=1:nl
    [q.rlnImageNo(i),q.rlnStackName{i}]=rlDecodeImageName(q.rlnImageName{i});
end;
disp(['Stack name: ' q.rlnStackName{1}]);
disp([num2str(nl) ' particles out of about ' num2str(q.rlnImageNo(end))]);
%%
% Count selected particles for each micrograph
snmi=numel(si.mi);
goodMiIndices=si.miIndex(q.rlnImageNo);
miGoodParticleNo=si.miParticle(q.rlnImageNo);
h=histogram(goodMiIndices,(0:snmi)+.5);
goodMiPicks=h.Values';

nPicks=0;
startEntry=1;
nmi=0;
miPicks=zeros(nEntries,1);
miDef=zeros(nEntries,1);
medianVar=zeros(nEntries,1);
medianGoodVar=zeros(nEntries,1);
medianAmp=zeros(nEntries,1);
medianGoodAmp=zeros(nEntries,1);
allVars=cell(nEntries,1);
allAmps=cell(nEntries,1);
allGoodVar=cell(nEntries,1);
allGoodAmp=cell(nEntries,1);

th0=2.88;  % base threshold
th1=1.4;  % quadratic term of threshold
mxv0=6;
mxv1=0.5;

okParticles=cell(nEntries,1);

for i=startEntry:nEntries
    if useAllMis
        mi=allMis{i};
        name=[infoDir mi.baseFilename 'mi.txt'];
    else
        name=['Info/' d(i).name];
        if ~strndcmp(name,'mi.txt')
            continue
        end;
        mi=ReadMiFile(name);
    end;
    nmi=nmi+1;
    d=mi.ctf(1).defocus;
    miDef(nmi)=d;
    thresh=th0+th1/d^2;
    maxVar=mxv0+mxv1*d;
    nMiParticles=size(mi.particle.picks,1);
    if nMiParticles>0
        flags=mi.particle.picks(:,3);
        amps=mi.particle.picks(:,5);
        vars=mi.particle.picks(:,8);
        sels=mi.particle.picks(:,10)>0; %%%%%%%%%%%%
        %             num=sum(flags>=16 & flags<48);
        okAmps=(flags==32 & amps>thresh & sels); % column vectors
        okVars=(flags==32 & vars<maxVar & sels);
        medianVar(nmi)=median(vars(okAmps));
        medianAmp(nmi)=median(amps(okVars));
        allVars{nmi}=vars(okAmps);
%  disp([num2str(nmi) ' ' num2str( max(allVars{nmi}))]);
        allAmps{nmi}=amps(okVars);
        
        rightMi=(goodMiIndices==nmi);
        miGoodParticles=miGoodParticleNo(rightMi);
        allGoodVar{nmi}=mi.particle.picks(miGoodParticles,8);
        medianGoodVar(nmi)=median(allGoodVar{nmi});
        allGoodAmp{nmi}=mi.particle.picks(miGoodParticles,5);
        medianGoodAmp(nmi)=median(allGoodAmp{nmi});
        
        okp=okAmps & vars<maxVar;
        num=sum(okp);
        okParticles{nmi}=okp;
        nPicks=nPicks+num;
        miPicks(nmi)=num;
    else
        num=0;
        okParticles{nmi}=[];
    end;
    if mod(nmi,stride)==0 || i>nEntries-3  % print out every 'stride' entries
        disp([num2str(nmi) ' ' name '  ' num2str([num nPicks])]);
    end;
end;
%%
figure(1);
miPicks=miPicks(1:nmi);
miDef=miDef(1:nmi);
figure(1);
subplot(413)
hist([miPicks goodMiPicks],100);
xlabel('Particles per micrograph');
ylabel('Frequency');

subplot(414);
hist(miDef,100);
xlabel('Defocus, um');
ylabel('Frequency');
subplot(412)
plot([miPicks goodMiPicks]);
ylabel('Particles per micrograph');
xlabel('Micrograph');
subplot(411);
plot(miDef);
ylabel('Defocus, \mum');
xlabel('Micrograph');

figure(2);
subplot(311);
plot(miDef,[miPicks goodMiPicks+.5],'.');
xlabel('Defocus, um');
ylabel('Number of picks');
axis([0 inf 0 30]);

subplot(312);
plot(miDef,[medianVar medianGoodVar],'.','markersize',10);
xlabel('Defocus, um');
ylabel('Median var');

subplot(313);
plot(miDef,[medianAmp medianGoodAmp],'.','markersize',10);
xlabel('Defocus, um');
ylabel('Median amp');


% figure(3);
% subplot(311);
% plot(miDef,[miPicks goodMiPicks+.5],'.');
% xlabel('Defocus, um');
% ylabel('Number of picks');
ds=1:.1:max(miDef);

varOffset=.1;
ampOffset=.02;
subplot(312);
cla;
hold on;
for i=1:nmi
     plot(repmat(miDef(i),numel(allVars{i}),1), allVars{i},'b.','markersize',5);
    plot(repmat(miDef(i),numel(allGoodVar{i}),1), allGoodVar{i}+varOffset,'r.','markersize',5);
end;
plot(ds,mxv0+mxv1*ds,'k-','linewidth',2);
hold off;

xlabel('Defocus, um');
ylabel('var');

subplot(313);
cla;
hold on;
for i=1:nmi
    plot(repmat(miDef(i),numel(allAmps{i}),1), allAmps{i},'b.','markersize',5);
    plot(repmat(miDef(i),numel(allGoodAmp{i}),1), allGoodAmp{i}+ampOffset,'r.','markersize',5);
end;
plot(ds,th0+th1*ds.^-2,'k-','linewidth',2);
hold off;

totalParticles=[sum(miPicks) sum(goodMiPicks)]
particlesPerImage=totalParticles/nmi 

return

