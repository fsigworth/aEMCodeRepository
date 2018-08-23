% rsCountParticles
% Count all the picks in a directory full of info files
maxEntries=inf;
stride=100;
d=dir('Info/');
nPicks=0;
nEntries=min(numel(d),maxEntries);
startEntry=5190;
nmi=0;
miPicks=zeros(nEntries,1);
miDef=zeros(nEntries,1);
for i=startEntry:nEntries
    
    name=['Info/' d(i).name];
    if strndcmp(name,'mi.txt')
        nmi=nmi+1;
        mi=ReadMiFile(name);
        miDef(nmi)=mi.ctf(1).defocus;
        if size(mi.particle.picks,1)>0
            flags=mi.particle.picks(:,3);
            num=sum(flags>=16 & flags<48);
            nPicks=nPicks+num;
            miPicks(nmi)=num;
        else
            num=0;
        end;
        if mod(nmi,stride)==0 || i>nEntries-3  % print out every 'stride' entries
            disp([num2str(nmi) ' ' d(i).name '  ' num2str([num nPicks])]);
        end;
    end;
end;
miPicks=miPicks(1:nmi);
miDef=miDef(1:nmi);
subplot(313)
hist(miPicks,100);
xlabel('Particles per micrograph');
ylabel('Frequency');
subplot(312)
plot(miPicks);
ylabel('Particles per micrograph');
xlabel('Micrograph');
subplot(311);
plot(miDef);
ylabel('Defocus, \mum');
xlabel('Micrograph');

