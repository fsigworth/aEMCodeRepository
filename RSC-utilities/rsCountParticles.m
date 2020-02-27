% rsCountParticles
% Count all the picks in a directory full of info files
useAllMis=1;
loadAllMis=0;

if loadAllMis
    MiLoadAll
end;

infoDir='Info/';
maxEntries=inf;
stride=1000;
if useAllMis
    load([infoDir 'allMis.mat']);
    nEntries=numel(allMis);
else
    d=dir(infoDir);
    nEntries=min(numel(d),maxEntries);
end;
%%
miNPicks=zeros(nEntries,1);
miDef=zeros(nEntries,1);
medianVar=zeros(nEntries,1);
nPicks=0;
startEntry=1;
nmi=0;


% for sq04_1
% Criteria for "good" particles
th0=2.88;  % base threshold
    th0=2.8; %%% new base.
th1=1.4;  % quadratic term of threshold
%    th1=2;  % quadratic term of threshold
mxv0=6;
    mxv0=7; %% var is not very sensitive
mxv1=0.4; % linear term of variance
maxFiltMax=3.8;
minDefocus=1.5;

maxNum=30;  % ignore any micrograph with more than this number of otherwise ok particles.

% % for sq07_2
% th0=3
% mxv0=7
% mxv1=0.7;
% 
% %for sq07_3
% mxv1=0.5

okParticles=cell(nEntries,1);

for i=startEntry:nEntries % loop over micrographs
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
    np0=size(mi.particle.picks,1);
    if np0>0
        flags=mi.particle.picks(:,3);
        amps=mi.particle.picks(:,5);
        vars=mi.particle.picks(:,8);
        filtMax=mi.particle.picks(:,11);
        defs=d+0*vars;
        %             num=sum(flags>=16 & flags<48);
        okAmps=(flags==32 & amps>thresh);
        okFiltMax=filtMax<maxFiltMax;
        okVars=vars<maxVar;
        okDefocus=defs>minDefocus;
        
        medianVar(nmi)=median(vars(okAmps));
        
        okp=okAmps & okFiltMax & okDefocus & okVars;
        okNum=sum(okp)<maxNum;
        okp=okp & okNum;
        num=sum(okp);

        disp([num2str(nmi) ' ' num2str([np0 num]) '  ' name]);
        
        okParticles{nmi}=okp;
        nPicks=nPicks+num;
        miNPicks(nmi)=num;
    else
        num=0;
        okParticles{nmi}=[];
    end;
%     if mod(nmi,stride)==0 || i>nEntries-3  % print out every 'stride' entries
%         disp([num2str(nmi) ' ' name '  ' num2str([num nPicks])]);
%     end;
end;
%
figure(2);
miNPicks=miNPicks(1:nmi);
miDef=miDef(1:nmi);
subplot(413)
hist(miNPicks,100);
xlabel('Particles per micrograph');
ylabel('Frequency');

subplot(414);
hist(miDef,100);
xlabel('Defocus, um');
ylabel('Frequency');
subplot(412)
plot(miNPicks);
ylabel('Particles per micrograph');
xlabel('Micrograph');
subplot(411);
plot(miDef);
ylabel('Defocus, \mum');
xlabel('Micrograph');


figure(3);
subplot(211);
plot(miDef,miNPicks,'o');
xlabel('Defocus, um');
ylabel('Number of picks');
subplot(212);
plot(miDef,medianVar,'o');
xlabel('Defocus, um');
ylabel('Median var');

totalParticles=sum(miNPicks)
particlesPerImage=totalParticles/nmi


return

%% Write modified info files.
doWrite=1;

allMis1=allMis;
allMis_old=allMis;
for i=1:nEntries
% for i=1:1
    mi=allMis{i};
    mi1=mi;
%     mi1.particle.picks=mi.particle.picks(okParticles{i},:);
%   add a 10th entry into picks, being a flag for valid particle.
    okp=single(okParticles{i});
        if numel(okp)>0
            mi1.particle.picks(:,10)=okp;
        end;
    allMis1{i}=mi1;
end;
if doWrite
    disp('saving allMis1.mat');
allMis=allMis1;
save([infoDir 'allMis1.mat'],'allMis');
allMis=allMis_old;  % restore the original one!
else
    disp('allMis1 created.');
end;
