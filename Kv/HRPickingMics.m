% HRPickMicrographs

if exist('isCluster','var') && isCluster
%     we have micInds already defined.
else
    micInds=8;
end;

% Locate particles on micrographs or in stacks
%

% % Fake data, start in this directory:
% cd ~/EMWork/Simulations/relion2/Fred2_HRPicking/
% starDir='Refine3D/job037/'
% starName='run_data.star'
% refDir=starDir;
% refName='run_class001.mrc'
% stackPrune=4;
% refCrop=96;
% bValue=0

% ----------
% W366F data start here:
% cd ~/EMWork/Yangyu/20210224_YW_sel/
cd '/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Yangyu/20210224_YW_sel'
% cd /Volumes/EMWork/Yangyu/20210224_YW_sel/
% cd ~/EMWork/20210224_YW_sel/

micDir='/Volumes/Drobo4/Yangyu_good_datasets/20210224/';
micDir=''

rootDir=AddSlash(pwd);
realData=1;
doCrossCorrelation=0;  %%%%%
ccUpsampling=1.5;
micHP=.002;
sharpHP=.06;

if realData
    starDir='Refine3D/job110/';
    starName='run_data.star';
else
    simDir='~/EMWork/20210224_YW_sel/Simulations/'
    cd(simDir);
    starDir='Stars1/'
    starName='particles_noVesPsi_varGroups_noOrigPix.star';
end;

% eigsName=[rootDir 'HRPicking/Eigs/Eigs_48.mat'];
eigsName=[rootDir 'HRPicking/Eigs/EigsTM_48.mat'];
useEigenimages=0;

projsName=[rootDir 'HRPicking/Eigs/projsTM.mat'];
projsName=[rootDir 'HRPicking/Eigs/projsComp_56.mat'];


outBasename=[rootDir 'HRPicking/CCsComp56_HP06'];
% To use the reconstructed map
% refDir='Postprocess/job171/';
% refName='postprocess.mrc';

bValue=0;
maxLines=1000;

skipStar=0;
if ~skipStar
    dataStarName=[starDir starName];
    disp(['Reading ' dataStarName ':']);
    [nms,dat]=ReadStarFile(dataStarName,1,maxLines);
    op=dat{1};
    d=dat{2};
    %%
    [micUNames,micUPtrs,micRPtrs]=unique(d.rlnMicrographName);
    nUMics=numel(micUNames); % number of mic names = no. of CTFs
    micUCounts=zeros(nUMics,1);
    for i=1:nUMics % we assume that micRPtrs are increasing stepwise
        micUCounts(i)=sum(micRPtrs==i);
    end;
end;
%% Collect some statistics
    nBins=100;
    micSigSampleMeds=zeros(nUMics,1);
    micSigSampleHists=zeros(nUMics,nBins);
    sigSampleBins=.5:nBins+.5;
    micMaxProbMeds=zeros(nUMics,1);
    micMaxProbHists=zeros(nUMics,nBins);
    maxProbBins=(0:nBins)/nBins;
    for j=1:nUMics
        micPtrs=micUPtrs(j):micUPtrs(j)+micUCounts(j)-1;
        q1=d.rlnNrOfSignificantSamples(micPtrs);
        micSigSampleMeds(j)=median(q1);
        micSigSampleHists(j,:)=histcounts(q1,sigSampleBins);
        q2=d.rlnMaxValueProbDistribution(micPtrs);
        micMaxProbMeds(j)=median(q2);
        micMaxProbHists(j,:)=histcounts(q2,maxProbBins);
    end;
    figure(9);
    subplot(221);
    bar(sigSampleBins(1:end-1)+.5,sum(micSigSampleHists,1));
    xlabel('No. significant samples');
    subplot(223);
    bar(maxProbBins(1:end-1)+.005,sum(micMaxProbHists,1));
    xlabel('Max probability');
    subplot(222);
    plot([micUCounts micSigSampleMeds]);
    legend('No. particles','Sign samples');
    subplot(224);
    plot(micMaxProbMeds);
    ylabel('Max prob medians');
    xlabel('Micrograph index');
    figure(8);
    plot3(micUCounts,micSigSampleMeds,1:nUMics,'o');
    xlabel('NumParticles');
    ylabel('NumSignificantSamples')
    zlabel('Index');
%     Some good inidices are 898, 1024, 1335
%     others (highest max prob) 1461 (25 parts), 1103 (88 parts)
% return
%%
% % Look at micrograph
% ds=3;
% iMic=16;
% while iMic>0
%     micLines=find(micRPtrs==iMic);
%     micLine=micLines(1);
%     micName=[micDir micUNames{iMic}];
%     disp(['Reading ' micName]);
%     [mic,s]=ReadMRC(micName);
% 
%     mOrig=RemoveOutliers(mic);
%     m0=GaussHP(mOrig,micHP*s.pixA);
%     sz1=NextNiceNumber(size(m0)/ds,11)*ds;
%     m0Mean=mean(m0(:));
%     m1=Crop(m0-mean(m0(:)),sz1);  % pad to a nice size.
% 
%     nd=round(size(m1)/ds);
%     md=Downsample(m1,nd);
%     mdf=GaussFilt(md,.2);
%     mdShift=(sz1-size(m0))/(2*ds);
%     imags(mdf);
%     txt=[num2str(iMic) ':  def ' num2str(d.rlnDefocusU(micLine)/1e4) '  nSignif  ' ...
%         num2str(micSigSampleMeds(iMic)) ...
%         '  maxProb ' num2str(micMaxProbMeds(iMic)) '  brt ' num2str(median(mOrig(:)))];
%     title(txt);
%     iMic=MyInput('Mic index? ',iMic+1);
% end;


%% Get a ctf for each micrograph
theMicInds=1335;
theMicInds=1103;
theMicInds=8;
% theMicInds=1461;
nUMics=numel(theMicInds);
activeMics=false(nUMics,1);
for i=1:nUMics
    iMic=theMicInds(i);
    micName=[micDir micUNames{iMic}];
    if exist(micName,'file')
        activeMics(i)=true;
    end;
end;
disp([num2str(sum(activeMics)) ' micrographs found.'])
activeMicInds=theMicInds(activeMics);
activeMicLines=micUPtrs(activeMicInds);

if numel(activeMicLines)<1
    disp('No micrographs found.');
    return
end;

disp('Making CTFs');
ct=rlStarLinesToCtf(nms,dat,activeMicLines);
nct=numel(ct);
for i=1:nct
    ct(i).B=bValue;
end;


%% Try finding particles on micrographs
if useEigenimages
    ei=load(eigsName);
    ds=ei.ds;
elseif doCrossCorrelation
    load(projsName);
end;
ds=3;

if sum(activeMics)<1
    disp('No micrographs found.')
    return
end;

figure(1);
for i=1:numel(activeMicInds)
    iMic=activeMicInds(i);
%     micName=micUNames{iMic};
    mLines=find(micRPtrs==iMic);
    micName=[micDir micUNames{iMic}];
    disp(['Reading ' micName]);
    [mic,s]=ReadMRC(micName);
    m0=RemoveOutliers(mic);
    m0=GaussHP(m0,micHP*s.pixA);

    sz1=NextNiceNumber(size(m0)/ds,11)*ds
    m0Mean=mean(m0(:));
    m1=Crop(m0-mean(m0(:)),sz1);  % pad to a nice size.

    nd=round(size(m1)/ds);
    md=Downsample(m1,nd);
    mdf=GaussFilt(md,.1);
    mdShift=(sz1-size(m0))/(2*ds);

    mdsh=SharpHP(m0,sharpHP*s.pixA*ds);
    fMic=fftn(-mdsh);
    ct(i).B=bValue;
    ct(i).pixA=s.pixA*ds;
    c=CTF(nd,ct(i));
    fMc=fftshift(fftn(md)).*c;
    mc=real(ifftn(ifftshift(fMc))); % ctf-modified, downsampled, HP filtered
    imags(mc);
    drawnow;
            outName=[outBasename num2str(i) '_' num2str(iMic) '.mat'];
    if doCrossCorrelation
        ncStep=1;
        disp('Cross correlations...')
        if useEigenimages
            ei.normCoeffs=ei.normCoeffs(1:ncStep:end,:);
            szr=size(ei.eigImgs,1);
            [mxVals,mxInds]=hrCCSearch(mc,ei,ccUpsampling);
        save(outName,'mxVals','mxInds','m1','mc','ct','c','s','ei');

        else
            szr=size(projs,1);
            [mxVals,mxInds]=hrProjSearch(mc,projs,ccUpsampling);
        save(outName,'mxVals','mxInds','m1','mc','ct','c','s');
        end;
        figure(2);
        imags(mxVals);
        return
    else
        load(outName);
    end;
end; % for i

return

%% %     cc2=GaussHP(mxVals,.005); %%%% in case we didn't do HP before.
dsv=2;
% Finding loop
dsv=2; %%%%%%
    cc2=mxVals;
    nPks=300;
    rMsk=30*dsv;
    nMsk=round(1.1*rMsk)*2;
    blankMask=1-fuzzymask(nMsk,2,rMsk);
    gix=zeros(nPks,1);
    giy=zeros(nPks,1);
    minVal=min(cc2(:));
    gpk=minVal*ones(nPks,1, 'single');
    for i=1:nPks
        [gpk(i),gix(i),giy(i)]=max2d(cc2);
        cc2=Mask(cc2,double([gix(i) giy(i)]),blankMask,(1-blankMask)*minVal);
%         if mod(i,25)==0
%             imags(cc2); drawnow;
%         end;
     end;
%     figure(3);
%     plot(gpk,'.-');
     gCoords=double([gix giy]);
%      save(outName,'mxVals','mxInds','gix','giy','gpk','m1','mc','ct','c','s');
%     disp(['Wrote ' outName]);
%
    %
% Draw boxes from the Relion coordinates
    rgnSize=100;
    gMax=60;  % no. of our CC peaks to show

    ctrReg=floor((rgnSize+1)/2);
    npl=numel(mLines);
    dsd=1;
    pX=d.rlnCoordinateX(mLines)+mdShift(1)*ds;
    pY=d.rlnCoordinateY(mLines)+mdShift(2)*ds;

    %  CC peak coordinates
    qX=gCoords(1:gMax,1)*dsv;
    qY=gCoords(1:gMax,2)*dsv;


% purge redundant particles
    dists=sqrt((pX-pX').^2+(pY-pY').^2);
    activeParts=true(npl,1);
    foundParts=true(npl,1);
    for i=1:npl
        for j=i+1:npl
            if dists(i,j)<rgnSize/2
                activeParts(j)=false;
            end
        end;
    end;
sum(activeParts)
micLines=mLines(activeParts);
np=sum(activeParts);

% Relion's purged, extracted coordinates
    pX=d.rlnCoordinateX(micLines)+mdShift(1)*ds;
    pY=d.rlnCoordinateY(micLines)+mdShift(2)*ds;
% Relion' particle shifts
    oShifts=d.rlnOriginXAngst(micLines)/op.rlnImagePixelSize(1);
    oShifts(:,2)=d.rlnOriginYAngst(micLines)/op.rlnImagePixelSize(1);

% Find true positive CC peaks
    ccDists=sqrt((pX-qX').^2+(pY-qY').^2);
%     truePos=false(np,1);
%     for j=1:np
%         if any(ccDists(j,:)<rgnSize,2)
%             truePos(j)=true;
%         end;
%     end;
truePos=any(ccDists<rgnSize,2);
nTruePos=sum(truePos)

    %     Make boxes for extracted particles, perhaps with Relion shifts
    [bX,bY,tX,tY]=MakeBoxDrawingVectors( ...
        ([pX pY]+0*oShifts)/dsd, ...
        op.rlnImageSize(1)/(2*dsd),0.8);
%     labels for the extracted particles
    tStrings=cell(np,1);
    for j=1:np
        tStrings{j}=num2str(micLines(j));
    end;

    %     Make boxes and labels for gMax CC peaks
    [gX,gY,wX,wY]=MakeBoxDrawingVectors( ...
        [qX qY]/dsd, op.rlnImageSize(1)/3*dsd,1);
    aStrings=cell(gMax,1);
     for j=1:gMax
%         if truePos(j)
%             ast='++';
%         else
            ast='';
%         end;
        aStrings{j}=[num2str(gpk(j),3) ast];
    end;

    figure(10); % --------------- Label the original image -------------
%    imags(GaussFilt(m1,.05));
%     imags(GaussFilt(mdsh,.1));
disLP=.1;
    imaga(imscale(GaussFilt(mdsh,disLP),256,.01));
    %         imags(mxVals.^2);
    hold on;
    plot(bX,bY,'color',[1 1 0],'linewidth',1); % extracted particles
    text(tX,tY,tStrings,'color',[1 1 0], ...
        'HorizontalAlignment','left',  'VerticalAlignment','top');
    plot(gX,gY,'color',[0.6 1 .4],'linewidth',1); % CC peaks
    text(wX,wY,aStrings,'color',[0.6 1 .4], ...
        'HorizontalAlignment','left',  'VerticalAlignment','top');
%         plot(gCoords(1:gMax,1)*dsv,gCoords(1:gMax,2)*dsv,'s','markersize',35,'color',[.8 1 .5]);
hold off;
    title([num2str(i) ':  ' num2str(np) ' original, ' num2str(gMax) ' CC peaks ' ...
    ' BP= [' num2str([sharpHP disLP]) ']']);
    drawnow;

    %     Get peak values near each particle
    peaks=zeros(np,1);
    coords=zeros(np,2);
    minDist=zeros(np,1);
    matchInd=zeros(np,1);
    for j=1:np
        pos=[pX(j) pY(j)]/dsv;
        rgn=ExtractImage(mxVals,round(pos), rgnSize);
        [peaks(j),xp,yp]=max2d(rgn);
        coords(j,:)=[xp yp]+round(pos)-ctrReg; % local peaks near extr. particle position
        matchDifs=sqrt(sum((repmat(coords(j,:),nPks,1)-gCoords).^2,2));
        [minDist(j),matchInd(j)]=min(matchDifs*dsv);
    end;
    goodDists=minDist<50;
    goodInds=matchInd(goodDists);

    %
    figure(11); % Show CC image
    imags(mxVals);
    dsv=2;
    hold on;
    plot(bX/dsv,bY/dsv,'color',[1 1 0]); % extr particles as boxes
    text(tX/dsv,tY/dsv,tStrings,'color',[1 1 0], ...
        'HorizontalAlignment','left',  'VerticalAlignment','top');

    %     plot matched local coords as circles
    plot(coords(:,1),coords(:,2),'yo','markersize',10); % circles: matched local coords

%     plot all CC peaks as squares
    plot(gCoords(1:gMax,1),gCoords(1:gMax,2),'s','markersize',20,'color',[.8 1 .4]);
%     plot(gCoords(goodInds,1),gCoords(goodInds,2),'s','markersize',20,'color',[1 .5 .5]);
    hold off;

    figure(12);
    plot(gpk,'.-');
    hold on;
    plot(goodInds,gpk(goodInds),'r+');
    hold off;


    %%
    % try to find particles

    %
    % ccs=zeros(s.mx,s.my,nRefs,'single');
    % for j=1:nRefs
    %         cpx=Crop(cProjs(:,:,j),[s.mx s.my]);
    %         ccs(:,:,j)=ifftshift(real(ifftn(fMic.*conj(fftn(cpx)))));
    % end;
    % cc=max(ccs,[],3);
    %
    %         figure(11);
    %         imags(cc);

    %
    %
    %
    %




