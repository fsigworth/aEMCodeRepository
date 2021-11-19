% HRPickingMicsFromMi.m
%  Reading mi files, compute cross correlations


cd( '/Volumes/D255/20181216/No5Graphene/sq05_1');

refDir='HRRefs/';
outDir='HRPicking/';
suffix='hp005';
realData=1;
doCrossCorrelation=1;  %%%%%
ds=3; % downsampling of refs
ccUpsampling=1.5;
% micHP=.002;
sharpHP=.005;
bValue=80;

projsName=[refDir 'projsComp_56.mat']
% % load(projsName);

miNames=f2FindInfoFiles;
maxMis=min(1,numel(miNames))

if doCrossCorrelation
    disp('Reading projections');
    load(projsName);
end;
for miInd=1:maxMis
    
    miName=miNames{miInd};
    disp(['Reading ' miName]);
    mi=ReadMiFile(miName);
    outName=[outDir mi.baseFilename '_' num2str(miInd) '_' suffix '.mat'];
    if doCrossCorrelation
        m1=meReadMergedImage(mi,0,'v');
        n1=size(m1);
        msk=hrpMaskDots(m1);
        m1h=SharpHP(m1,sharpHP*mi.pixA);
        sd=.1*std(m1h(:));
        m1h=m1h.*(1-msk)+sd*msk.*randn(n1);

        nd=round(n1/ds); % downsampled to CC size
        
        mdh=Downsample(m1h,nd);
        imags(mdh);
        
        m1u=meReadMergedImage(mi,0); % unsubtracted image

        
        pixAd=mi.pixA*ds;
        
        cPars=mi.ctf(1);
        cPars.pixA=mi.pixA*ds;
        cPars.B=bValue;
        c=CTF(nd,cPars);
        
        %% Try finding particles on micrographs        
        fMc=fftshift(fftn(mdh)).*c;
        mc=real(ifftn(ifftshift(fMc))); % ctf-modified, downsampled, HP filtered
        
        figure(1); clf;
        imags(mc);
        drawnow;
        nProjs=size(projs,3);
        %         crProjs=projs(:,:,1:100);
        %                 nProjs=size(crProjs,3);
        dsu=ds*ccUpsampling;
        
            disp(['Cross correlation with ' num2str(nProjs) ' references...'])
            
            [mxVals,mxInds,sumCCs,sumSquares]=hrProjSearch(mc,projs,ccUpsampling);
            
            corrMxVals=(mxVals-sumCCs/nProjs)./sqrt(sumSquares/nProjs-(sumCCs/nProjs).^2);
            CheckAndMakeDir(outDir);
            save(outName,'mxVals','mxInds','sumCCs','sumSquares',...
                'corrMxVals','m1','m1u','mc','mi','dsu');
    else
        disp(['Loading ' outName]);
        load(outName);
  
    end;
%%
        figure(1); clf;
        mysubplot(131);
        imags(m1u);
           mysubplot(132);
        imags(mxVals);
        mysubplot(133);
        imags(corrMxVals);
        drawnow;
end; % for miInd
return

%% %     cc2=GaussHP(mxVals,.005); %%%% in case we didn't do HP before.
% Finding loop
dsv=ds*ccUpsampling;
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
    gMax=10;  % no. of our CC peaks to show

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
    truePos=false(np,1);
    for j=1:np
        if any(ccDists(j,:)<rgnSize,2)
            truePos(j)=true;
        end;
    end;
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
        if truePos(j)
            ast='++';
        else
            ast='';
        end;
        aStrings{j}=[num2str(gpk(j),3) ast];
    end;

    figure(10); % Label the original image
    imags(GaussFilt(m1,.05));
    %         imags(mxVals.^2);
    hold on;
    plot(bX,bY,'color',[1 1 0],'linewidth',1.5); % extracted particles
    text(tX,tY,tStrings,'color',[1 1 0], ...
        'HorizontalAlignment','left',  'VerticalAlignment','top');
    plot(gX,gY,'color',[0.6 1 .4],'linewidth',1); % CC peaks
    text(wX,wY,aStrings,'color',[0.6 1 .4], ...
        'HorizontalAlignment','left',  'VerticalAlignment','top');
%         plot(gCoords(1:gMax,1)*dsv,gCoords(1:gMax,2)*dsv,'s','markersize',35,'color',[.8 1 .5]);
hold off;
    title([num2str(i) ':  ' num2str(np) ' original, ' num2str(gMax) ' CC peaks']);
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




