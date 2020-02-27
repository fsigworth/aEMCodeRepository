% rlShowParticlesOnMicrographs
% Based on rsCountPicks2
% Read a (possibly merged) star file. Read in one or more corresponding
% si files. For a given micrograph, show the particles based on classes.

return % don't allow running from the beginning....

% % Run this in the Stack directory.
cd('/Volumes/Extreme1T/20181218/No5/sqall/Stack3cp');
load('../../../workDirs.mat');

localWorkDirs=workDirs;
disp('Working directories:');
for i=1:numel(workDirs)
    localWorkDirs{i}=workDirs{i}(26:end);
    disp(localWorkDirs{i});
end;


% starName='Dec20_02all.star';
% starName='Class2D_job008_data.star';
starName='Class3D034/run_it025_data.star';
%         rlnGroups         our tag
% dotClasses=[7 28 31 43 45]; % 5 : dots
% edgeDotClasses=[18 19 20 24 0]; % 4 : edge dots
% mbnClasses=[13 0 0 0 0];   % 3 : membranes
% partClasses=[11 50 0 0 0]; % 2 : particles
% junkClasses=[];            % 1 : noise
% sieve=[partClasses; mbnClasses; edgeDotClasses; dotClasses];
% assign a tag number by ca=find(any(sieve==x,2)) + 1
% if numel(ca)<1
%     ca=1;
% end;
% TagAssignment=ca;

% Get the full particle set
disp(['Reading ' starName]);
[~,dat]=ReadStarFile(starName);
d=dat{1};

% % Get the 3D Class particles, a subset
% disp(['Reading ' starName1]);
% [~,dat1]=ReadStarFile(starName1);
% d1=dat1{1};
% 
% %%
% % Find pointers from star1 to star
% iName=string(d.rlnImageName); % string arrays run twice as fast.
% iName1=string(d1.rlnImageName);
% disp('Matching star lines');
% ptrs=zeros(n1,1);
% 
% for i=1:n1 % takes about 16 ns per comparison...
%     ptrs(i)=find(strcmp(iName1(i),iName),1);
% end;

% Extract the si file names
disp('Extracting si names');
nl=numel(d.rlnImageName);
mName=cell(nl,1);
pIndex=zeros(nl,1);
for i=1:nl
    [pIndex(i),mName{i}]=rlDecodeImageName(d.rlnImageName{i});
end;
%% Load the si structures into sis{}
[micNames,micInds,indsAll]=unique(mName); % micInds is the starting line of each si set.
nsi=numel(micNames);
siNames=micNames;
sis=cell(nsi,1);
for i=1:nsi
    siNames{i}=[micNames{i}(1:end-9) 'i.mat'];
    disp(['Reading ' siNames{i}]);
%     ok=exist(siNames{i},'file');
    s=load(siNames{i});
    sis{i}=s.si;
end;
%%
% Assemble the information about each particle
pi=struct; % particle info
pi.defocus=(d.rlnDefocusU+d.rlnDefocusV)/2e4; % mean defocus in um
pi.coords=zeros(nl,14,'single');
pi.pickerPars=zeros(nl,4,'single');
pi.miIndex=zeros(nl,1,'single');
pi.label=zeros(nl,2,'single');
pi.baseName=cell(nl,1);
pi.siIndex=indsAll;

parsText={'amp'; 'var'; 'filtMax'; 'filtMaxOfResidual'};
for i=1:nl
% pick up the relevant mi structure and parameters
    j=pIndex(i); % index into stack
    si=sis{pi.siIndex(i)}; % index into sis{}
    pi.miIndex(i)=si.miIndex(j); % index into si.miIndex
    mi=si.mi{si.miIndex(j)};     % mi
    pi.miParticle(i)=si.miParticle(j); % miParticle
    pi.baseName{i,:}=mi.baseFilename;
    coords=mi.particle.picks(pi.miParticle(i),:);
    nc=numel(coords);
    pi.coords(i,1:nc)=coords(1:nc);
    pi.pickerPars(i,:)=coords(1,[9 8 11 12]).*[1e4 1 1 1];
    %     pickerPars are amp, var, filtMax, filtMaxOfResidual
% %     x=d.rlnClassNumber(i);
% %     ca=find(any(sieve==x,2)) + 1;
% % %     labels are: 
% % 
% % %     1: noise; 2: particles; 3: membranes; 4: edge dots; 5: strong dots.
% % if numel(ca)<1
% %     pi.label(i,1)=1; % no assignment says it's noise.
% % end;
% % %     pi.label(i,1)=min(4,ca); % combine classes 4 and 5 (both dots)
% %     pi.label(i,1)=ca;
end;
% classLabels=["noise" "particles" "membranes" "edge dots" "dots" "best3D"];

% Set the labels for the 3D classes
% pi.label(ptrs,2)=d1.rlnClassNumber;

%%
fc=.05; % A^-1
markerSize=40;
partInd=1;
colors=[repmat([1 .5 0],5,1); .1 1 .3 ];
nColors=size(colors,1);
startPath='../../../../'; 
% Get a micrograph

showSubtracted=0;

%% partInd=myinput('particle index',partInd);
while partInd>0
    
    if partInd <= numel(pi.siIndex)
        
        % Pick the image
        siInd=pi.siIndex(partInd);
        miInd=pi.miIndex(partInd);
        mi=sis{siInd}.mi{miInd};
        baseImageName=mi.baseFilename;
        pixA=mi.pixA*mi.imageSize(1)/size(mi,1);
        workPath=localWorkDirs{siInd};
        imgName=[startPath workPath 'Merged/' baseImageName 'ms.mrc'];
        imgNamev=[startPath workPath 'Merged/' baseImageName 'mvs.mrc'];
        if showSubtracted
            imgName=imgNamev;
        end;
        disp(imgName);
        [m,s]=ReadMRC(imgName);
        imaga(imscale(GaussFilt(m,fc*s.pixA),256,1e-4));
        
        imageParticles=find(pi.siIndex==siInd & pi.miIndex==miInd);  % find all the particles
        nParticles=numel(imageParticles);
        for i=1:nParticles
            px=pi.coords(imageParticles,1)/4+1;
            py=pi.coords(imageParticles,2)/4+1;
            class=d.rlnClassNumber(imageParticles);
        end;
        hold on
        for i=1:nColors
            ok=class==i;
            if sum(ok)>0
                plot(px(ok),py(ok),'s','markerSize',markerSize,'color',colors(i,:));
            end;
        end
        hold off
        title([num2str(partInd) '    ' workPath '     ' baseImageName '  ' ...
            num2str(mi.ctf.defocus,3) ' um'],'interpreter','none');
        oldPartInd=partInd;
        partInd=partInd+max(1,nParticles);
        
    else
        disp('Out of range.');
        partInd=numel(pi.siIndex);
    end;
    
    partInd=MyInput('Particle index? ',partInd);
    if partInd<0
        partInd=oldPartInd;
        showSubtracted=1;
    else
        showSubtracted=0;
    end;
    
end;
disp('done.');

return


%% Get vesicle size distribution fpr vesicles containing good particles
radii=zeros(nl,1);
amps=zeros(nl,1);
for i=1:nl
    mi=sis{pi.siIndex(i)}.mi{pi.miIndex(i)};
    radii(i)=mi.vesicle.r(pi.coords(i,4),1)*mi.pixA;
    amps(i)=mi.vesicle.s(pi.coords(i,4),1);
end;

oks=amps>.0004 & amps <.002 & radii>50 & radii <250;
figure(3);
subplot(211);
hist(2*radii(oks),200);
xlabel('Vesicle diameter, A');
ylabel('Frequency');
subplot(212);
hist(amps(oks),200);
xlabel('Vesicle contrast, radians/2VÅ');
ylabel('Frequency');

%% Get the global size distribution
allRadii=[];
allAmps=[];
for i=1:nsi
    si=sis{i};
    for j=1:numel(si.mi);
        mi=si.mi{j};
        nv=numel(mi.vesicle.x);
        if nv>0
        allRadii(end+1:end+nv,1)=mi.vesicle.r(:,1)*mi.pixA;
        allAmps(end+1:end+nv,1)=mi.vesicle.s(:,1);
        end;
    end;
end;            
figure(4);
okRadii=allAmps>.0008 & allRadii>50 & allRadii <300 
subplot(211)
hist(2*allRadii(okRadii),300);
okAmpss=allAmps>0 & allAmps < .002;
subplot(212);
hist(allAmps(okAmps),300);
%%

workPath=localWorkDirs{pi.siIndex(partInd)};
imgName=[startPath workPath 'Merged/' sis{pi.siIndex(partInd)}.mi{pi.miIndex(partInd)}.baseFilename 'ms.mrc'];
disp(imgName);
m=ReadMRC(imgName);
pos=pi.coords(partInd,1:2)/4+1;

imags(m);
hold on;
plot(pos(1),pos(2),'gs','markerSize',30);
hold off;

return
%%


micInds=1:4;
mul=[1 .5 1 1]; % scaling of picker pars
figure(2);
xs=0.5:.02:8;
for j=1:4 % column: picker parameter
    for i=1:6 % row: class label
        subplot(6,4,j+4*(i-1));
        if i<6
            oks=pi.label(:,1)==i;
        else
            oks=pi.label(:,2)==i-5;
        end;
        %     disp(sum(oks));
        oksAll=oks & any(indsAll==micInds,2);
        oksConstr=oksAll & pi.pickerPars(:,3)<3.8 & pi.pickerPars(:,1)<5 ...
            & pi.pickerPars(:,1)>2.8 & pi.defocus>1.5;
        hist(pi.pickerPars(oksConstr,j)*mul(j),xs);
        if i==1
            title([num2str(j) ' ' parsText{j}]);
        end;
        if j==1;
            ylabel(classLabels(i));
            text(max(xs)*0.8,0,num2str([sum(oksConstr) sum(oksAll)]),'VerticalAlignment','bottom');
            if i==6
                xlabel(micNames(micInds),'interpreter','none');
            end;
        end;
    end;
end;










return;

%----------------------------------------------------


% Count all the picks in a directory full of info files. Compare these with
% particles in a star file from a selection operation.





useAllMis=1;
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
nl=numel(q.rlnImageName); % number of lines in star file
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

