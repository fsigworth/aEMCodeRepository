% HRPickMicrographs

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
rootDir=AddSlash(pwd);
realData=1;
doCrossCorrelation=0;
upsampling=1.5;

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
projsName=[rootDir 'HRPicking/Eigs/projsTM.mat'];
useEigenimages=0;


outName=[rootDir 'HRPicking/ETM_B0_W366F_CCs.mat'];
% To use the reconstructed map
% refDir='Postprocess/job171/';
% refName='postprocess.mrc';

bValue=00;
maxLines=1000;

skipStar=0;
if ~skipStar
    dataStarName=[starDir starName];
    disp(['Reading ' dataStarName '...']);
    [nms,dat]=ReadStarFile(dataStarName,1,maxLines);
    op=dat{1};
    d=dat{2};
    %%
    [micUNames,micUPtrs,micRPtrs]=unique(d.rlnMicrographName);
    nUMics=numel(micUNames); % number of mic names = no. of CTFs
end;
%% Get a ctf for each micrograph
activeMics=false(nUMics,1);
for i=1:nUMics
    if exist(micUNames{i},'file')
        activeMics(i)=true;
    end;
end;
disp([num2str(sum(activeMics)) ' micrographs found.'])
activeMicLines=micUPtrs(activeMics);

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
activeMicInds=find(activeMics);
figure(1);
for i=1
    iMic=activeMicInds(i);
    micName=micUNames{iMic};
    micLines=find(micRPtrs==iMic);
    % micName='Micrographs1/mic000.mrc'
    [mic,s]=ReadMRC(micName);
    m0=RemoveOutliers(mic);
    sz1=NextNiceNumber(size(m0)/ds,11)*ds
    m0Mean=mean(m0(:));
    m1=Crop(m0-mean(m0(:)),sz1);  % pad to a nice size.


    %     micParticles=find(micRPtrs==i);
    %     np=numel(micParticles);
    %     if np<1
    %         continue;
    %     end;
    nd=round(size(m1)/ds);
    md=Downsample(m1,nd);
    mdf=GaussFilt(md,.1);
    %     imags(mdf);
    fMic=fftn(-mic);
    ct(i).B=bValue;
    ct(i).pixA=s.pixA*ds;
    c=CTF(nd,ct(i));
    fMc=fftshift(fftn(md)).*c;
    mc=real(ifftn(ifftshift(fMc)));
    imags(mc);
    drawnow;
    %     imags(c);
    if doCrossCorrelation
        ncStep=1;
            disp('Cross correlations...')
        if useEigenimages
            ei.normCoeffs=ei.normCoeffs(1:ncStep:end,:);
            [mxVals,mxInds]=hrCCSearch(mc,ei,upsampling);
        else
            [mxVals,mxInds]=hrProjSearch(mc,projs,upsampling);
        end;
        figure(2);
        imags(mxVals);

        % Find the highest peaks
        %%
        cc2=mxVals;

        nPks=600;
        rMsk=25;
        nMsk=round(1.1*rMsk)*2;
        blankMask=1-fuzzymask(nMsk,2,rMsk);
        gix=zeros(nPks,1,'int32');
        giy=zeros(nPks,1,'int32');
        minVal=min(cc2(:));
        gpk=minVal*ones(nPks,1, 'single');
        for i=1:nPks
            [gpk(i),gix(i),giy(i)]=max2d(cc2);
            cc2=Mask(cc2,double([gix(i) giy(i)]),blankMask,(1-blankMask)*minVal);
            if mod(i,25)==0
                imags(cc2); drawnow;
            end;
        end;
        figure(3);
        plot(gpk);
        gCoords=double([gix giy]);
        save(outName,'mxVals','mxInds','gix','giy','gpk','m1','mc','ct','c','s');
        disp(['Wrote ' outName]);
    else
        load(outName);
    end;
    return

    %%
    rgnSize=100;
    ctrReg=floor((rgnSize+1)/2);
    np=numel(micLines);
dsd=1;
        pX=d.rlnCoordinateX(micLines);
        pY=d.rlnCoordinateY(micLines);
        [bX,bY,tX,tY]=MakeBoxDrawingVectors( ...
            [pX pY]/dsd, ...
            op.rlnImageSize(1)/(2*dsd),0.8);
        tStrings=cell(np,1);
        for j=1:np
            tStrings{j}=num2str(micLines(j));
        end;
        figure(10);
        imags(GaussFilt(m1,.1));
%         imags(mxVals.^2);
        hold on;
        plot(bX,bY,'color',[1 1 0]);
        text(tX,tY,tStrings,'color',[1 1 0], ...
            'HorizontalAlignment','left',  'VerticalAlignment','top');
        hold off;
        title(i);
        drawnow;
%
        figure(11);
        imags(mxVals);
        dsv=2;
        hold on;
        plot(bX/dsv,bY/dsv,'color',[1 1 0]);
        text(tX/dsv,tY/dsv,tStrings,'color',[1 1 0], ...
            'HorizontalAlignment','left',  'VerticalAlignment','top');
        hold off;
%     Get peak values near each particle    
peaks=zeros(np,1);
coords=zeros(np,2);
minDist=zeros(np,1);
matchInd=zeros(np,1);
for j=1:np
    pos=[pX(j) pY(j)]/dsv;
    rgn=ExtractImage(mxVals,round(pos), rgnSize);
    [peaks(j),xp,yp]=max2d(rgn);
    coords(j,:)=[xp yp]+round(pos)-ctrReg;
    matchDifs=sqrt(sum((repmat(coords(j,:),nPks,1)-gCoords).^2,2));
    [minDist(j),matchInd(j)]=min(matchDifs);
end;
goodDists=minDist<20;
goodDists=true(nPks,1);
hold on;
plot(coords(:,1),coords(:,2),'yo','markersize',10);
plot(gCoords(goodDists,1),gCoords(goodDists,2),'s','markersize',20,'color',[1 .5 .5]);
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
end;





