% rlCheckDefocus
% Get defocus and other stats for micrographs in a dataset
% % Put up a file selector for *si.mat,
% [siName, siPath]=uigetfile('*si.mat','Select si file');
% if isnumeric(stPath)  % user has clicked Cancel
%     return
% else
%     cd(siPath);
% end;
%
doLoadFiles=1;
doFileSelector=0;

if ~exist('doLoadFiles','var') || doLoadFiles
    if doFileSelector
    disp('Getting a data.star file from the classification.');
    [stName, stPath]=uigetfile('*.star','Select data star file');
    if isnumeric(stPath)  % user has clicked Cancel
        return
    else
        cd(stPath);
    end
    else
        
    stPath='/Users/fred/EMWork/Hideki/161101/KvLipo122_4b/Stack2/Class3D/job007';
    cd(stPath);
    stName='run_it025_data.star';
    end;
    [ba1,pa1]=ParsePath(stPath);
    [ba2,pa2]=ParsePath(ba1);
    [ba3,pa3]=ParsePath(ba2);
    localPath=[pa3 pa2 pa1];
    %%
    disp(['Reading ' stName '...']);
    [stBlocks,stData,ok]=ReadStarFile(stName);
    disp('done.');
    %%
    dat=stData{1};
    nim=numel(dat.rlnDefocusU);
    if isfield(dat,'rlnClassNumber')
        cls=dat.rlnClassNumber;
    else
        cls=zeros(nim,1);
    end
    def=(dat.rlnDefocusU + dat.rlnDefocusV)/2e4;  % defocus in um
    pInds=zeros(nim,1);
    for i=1:nim
        pInds(i)=rlDecodeImageName(dat.rlnImageName{i});
    end
%     pInds(i) is the original si file index for the ith selected particle.
    defRange=.5:.1:5.5;
    h0=hist(def,defRange);
    % figure(2);
    plot(defRange,h0);
    nTotal=nim
    % clsSelections=[8 9 25 42 55 56 58 66 68 73 94];  % nicest of 100.
    % clsSelections=90;
    nCls=max(cls);
    % sel=false(nim,1);
    % for i=1:nim
    %     sel(i)=any(clsSelections==cls(i));
    % end;
    
    % Make a histogram for each class
    h=zeros(numel(defRange),nCls);
    for i=1:nCls
        h(:,i)=hist(def(cls==i),defRange);
        h(:,i)=h(:,i)/sum(h(:,i));  % normalization
    end;
    
    % if max(h1)>0
    % hScale=max(h1)/max(h0);
    % else
    %     hScale=1;
    % end;
    % plot(defRange,[h0'*hScale h1']);
    figure(1);
    plot(defRange,h);

        %%
    
        if doFileSelector
            
            disp('Getting a si.mat file from the original images');
            [siName, siPath]=uigetfile('*si.mat','Select stack info file');
            if isnumeric(stPath)  % user has clicked Cancel
                return
            end;
            
        else
            siPath='/Users/fred/EMWork/Hideki/161101/KvLipo122_4b/Stack2/';
            siName='sq81_1p192m2+tsi.mat';
        end;
        
    load([AddSlash(siPath) siName]);
%%
    
    % Make a timecourse for each class
    nInd=max(pInds);  % go back to original image number.
    tCls=zeros(nInd,nCls);
    pCls=zeros(nim,nCls);
    for i=1:nCls
        pCls(:,i)=dat.rlnClassNumber==i;
        tCls(pInds,i)=cls==i;
    end;
    %
    points=[1 12380 16680 18890 22600 27610 30740 45000 63000 nim-3000];  % interesting particle numbers
    figure(2);
    subplot(211);
%     plot(GaussFilt(tCls,.002,1));
    gCls=GaussFilt(pCls,.002,1);
    plot(gCls);
    hold on;
    plot(points,gCls(points,2),'k+','markersize',10,'linewidth',1.5);
    
    for i=1:numel(points)
        txt=si.mi{si.miIndex(pInds(points(i)))}.baseFilename;
        text(points(i),-.1,{txt(1:11),txt(13:end)},'rotation',90,'interpreter','none','verticalalignment','top');
    end;
    hold off;
    axis([0 inf 0 inf]);
    grid on;
    xlabel('Particle number');
    ylabel('Class probability');
    title(localPath,'interpreter','none');
    LegendAuto(nCls);
 
    subplot(212)
%     gDef=GaussFilt(def,.002);
    gDef=def;
    plot(gDef);
    hold on;
    plot(points,gDef(points),'k+','markersize',10);
    hold off;
        xlabel('Particle number');
    ylabel('Defocus, um');
    grid on
    axis([0 inf 0 inf]);
end;

    
%% get micrograph index for each selected particle
disp('Indexing micrographs');
nmi=numel(si.mi);
mInd=zeros(nim,1);
cParticles=cell(nmi,1);
totalNmp=0;
miInds=si.miIndex(pInds);
for i=1:nmi
%     Search for a match between mi.baseName and all the particles'
%     micrograph names.
%     bName=si.mi{i}.baseFilename;
%     q=strncmp(dat.rlnMicrographName,bName,numel(bName));
    q=miInds==i;
    nmp=sum(q);  % number of micrograph particles
    if nmp>0
        mInd(q)=i;
        cParticles{i}=find(q);
        totalNmp=totalNmp+nmp;
    end;
end;


%% compute statistics for each micrograph

mCls=zeros(nmi,nCls);  % Class distribution for micrograph
mDef=zeros(nmi,1);
mAmp=zeros(nmi,1);     % mean amplitude
mClsAAmp=zeros(nmi,nCls);  % mean class abs amplitude for each micrograph
mShift=zeros(nmi,1);   % mean shift magnitude
mnParts=zeros(nmi,1);  % number of particles from micrograph
tnParts=zeros(nmi,1);
for i=1:nmi  % For each micrograph
    cParts=cParticles{i}; % the relevant particles in star file
    nmp=numel(cParts);
    mnParts(i)=nmp;  % micrograph no. particles
    for j=1:nCls
        mCls(i,j)=sum(cls(cParts)==j)/max(1,nmp);  % class distribution
    end;
    sParts=pInds(cParts); % get si particle indices
    mParticles=si.miParticle(sParts);
    mi=si.mi{i};
    if isfield(mi.particle,'picks')
        tnParts(i)=size(mi.particle.picks,1);
    end;
    mDef(i)=mi.ctf(1).defocus;
    if nmp>0
        if abs(mDef(i)-def(cParts(1)))>1e-3  % check that this agrees with si file
            disp([mDef(i) def(cParts(1))]);
        end;
        mShift(i)=sqrt(sum(mi.frameShifts{1}(1,:).^2));
        cpInds=si.miParticle(sParts);  % particles in the ith micrograph
        vesAmps=si.sVesicle(sParts);
vesicleInds=mi.particle.picks(cpInds,4);
okVes=(vesicleInds>0);
vesicleAmps=mi.vesicle.s(vesicleInds(okVes),1,1);
        particleVAmps=mi.particle.picks(cpInds(okVes),5);
        mAmp(i)=mean(particleVAmps.*vesicleAmps);
        pks=mi.particle.picks(cpInds,:);  % particle picks row
        if size(pks,2)>8 && any(pks(:,9)>0)
            mAmp(i)=mean(pks(:,9));
            for j=1:nCls
                mClsAAmp(i,j)=sum((cls(cParts)==j).*pks(:,9))/max(1,sum(cls(cParts)==j));
            end;
        end;
     end;
end;
%% show statistics per micrograph
mPoints=double(si.miIndex(pInds(points)));

figure(4);
nr=4;
subplot(nr,1,1);
fmCls=GaussFilt(mCls,.1,1);
% plot(mCls,'-');
plot(fmCls,'-');
hold on;
plot(mPoints,fmCls(mPoints,2),'k+','markersize',10,'linewidth',1.5);
for i=1:numel(mPoints)
    txt=si.mi{mPoints(i)}.baseFilename;
    text(mPoints(i),0.6,{txt(1:11),txt(13:end)},'rotation',90,'interpreter','none','verticalalignment','top');
end;
hold off;

% plot(GaussFilt(mCls,.05,1));
xlabel([localPath '   Micrograph'],'interpreter', 'none');
ylabel(['Class in ' stName],'interpreter','none')
LegendAuto(nCls);
% title(localPath,'interpreter','none');

subplot(nr,1,2)
plot(mAmp);
xlabel('Micrograph');
ylabel('Particle amplitude');

subplot(nr,1,3);
plot(mDef);
xlabel('Micrograph');
ylabel('Defocus 1');

subplot(nr,1,4);
% plot(GaussFilt(mShift,.05));
plot(mShift);
xlabel('Micrograph');
ylabel('Initial shift, pixels');


cls1=3;
cls2=2;



figure(5);  % scatterplots
cl4=[mCls(:,cls1)/mean(mCls(:,cls1)) mCls(:,cls2)/mean(mCls(:,cls2))];
% cl4=mCls(:,[cls1 cls2]);
% cl4(:,[ ])=0;
ylbl=['Class ' num2str(cls1) ' (blue) ' num2str(cls2) ' (red)'];

subplot(2,2,1);
plot(GaussFilt(mShift,.05),cl4,'.');
xlabel('Initial shift');
ylabel(ylbl);

subplot(2,2,2);
plot(mAmp.^2,cl4,'.');
xlabel('particle power');
ylabel(ylbl);

subplot(2,2,3);
plot(mDef,cl4,'.');
xlabel('Defocus');
ylabel(ylbl);


return

%% Compute particle amplitudes
startIndex=1224;
endIndex=startIndex-1+500;
basePath='/Users/fred/EMWork/Hideki/160909/KvLipo121_2retrack';
for i=startIndex:endIndex
    si.mi{i}=rspInsertAbsAmplitudes(si.mi{i},basePath);
end;



