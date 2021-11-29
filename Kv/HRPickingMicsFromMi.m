% HRPickingMicsFromMi.m
%  Reading mi files, compute cross correlations


cd( '/Volumes/D255/20181216/No5Graphene/sq05_1');

refDir='HRRefs/';
outDir='HRPicking/';
suffix='pw_B80_AlphaWeak';
realData=1;

doLoadProjections=1;
doCrossCorrelation=1;  %%%%%
finalFigIndex=1;

ds=3; % downsampling of refs
ccUpsampling=1.5;
% micHP=.002;
sharpHP=.005;
bValue=120;

% projsName=[refDir 'projsComp_56.mat']
% projsName=[refDir 'projsComp__2.5_2.5_4_56.mat']
% projsName=[refDir 'projsTM__4_4_4_48.mat']
projsName=[refDir 'projsCompWeak__4_4_4_56.mat']

miNames=f2FindInfoFiles;
maxMis=min(1,numel(miNames))

if doLoadProjections && doCrossCorrelation
    disp('Reading projections');
    load(projsName);
end;

maxMis=1
for miInd=1:maxMis

    miName=miNames{miInd};
    disp(['Reading ' miName]);
    mi=ReadMiFile(miName);
    outName=[outDir mi.baseFilename '_' num2str(miInd) '_' suffix '.mat'];
    if doCrossCorrelation
        m1=meReadMergedImage(mi,0,'v');
        m1u=meReadMergedImage(mi,0); % unsubtracted image
        n1=size(m1);

        % Create the dot mask
        %         msk=hrpMaskDots(m1);
        mf1=GaussFilt(GaussHP(m1,.003),.05);
        msk1=GaussFilt(mf1<-.11, .05)>.02;
        msk=(GaussFilt(msk1,.02)>.01);

        m1h=SharpHP(m1,sharpHP*mi.pixA);
        sd=.1*std(m1h(:));
        m1h=m1h.*(1-msk)+sd*msk.*randn(n1); % add random numbers outside the mask

        figure(2);
        pwf1=hrpMakePWFilter(m1,mi.pixA,1);
        m1pw=real(ifftn(fftn(m1h).*ifftshift(pwf1)));

        nd=round(n1/ds); % downsampled to CC size

        mdh=Downsample(m1pw,nd);
        pixAd=mi.pixA*ds;

        figure(3);
        mysubplot(221);
        imags(GaussFilt(m1u,.1*mi.pixA));

        mysubplot(222);
        imags(GaussFilt(mdh,.1*pixAd));
        title('Pre-whitened');

        cPars=mi.ctf(1);
        cPars.pixA=mi.pixA*ds;
        cPars.B=bValue;
        c=CTF(nd,cPars);

        %% Try cross-correlating
        fMc=fftshift(fftn(mdh)).*c;        
        mc=real(ifftn(ifftshift(fMc))); % ctf-modified, downsampled, HP filtered

        mysubplot(223);
        imags(mc);
                title('CTF-multiplied');
        drawnow;
        nProjs=size(projs,3);
        %         crProjs=projs(:,:,1:100);
        %                 nProjs=size(crProjs,3);
        dsu=ds*ccUpsampling;
        if MyInput('Ready to do cross correlations',0) ==1
            disp(['Cross correlation with ' num2str(nProjs) ' references...'])

            [mxVals,mxInds,sumCCs,sumSquares]=hrProjSearch(mc,projs,ccUpsampling);

            corrMxVals=(mxVals-sumCCs/nProjs)./sqrt(sumSquares/nProjs-(sumCCs/nProjs).^2);
            CheckAndMakeDir(outDir);
            save(outName,'mxVals','mxInds','sumCCs','sumSquares',...
                'corrMxVals','m1','m1u','mc','mi','dsu','projsName','outName');
        else
            disp('Stopped.');
        end
    else
        disp(['Loading ' outName]);
        load(outName);

    end;

%%

    n1=size(m1);
    figure(1); clf;
    mysubplot(221);
    imags(m1u);
    mysubplot(222);
    imags(m1);
    mysubplot(223);
    imags(mxVals);
    mysubplot(224);
    imags(corrMxVals);
    drawnow;

%


figure(5+finalFigIndex); clf;
%%

bkg=rot90(imscale(Downsample(GaussFilt(m1u,.1),1920),128));
imaga(bkg)
rgb=repmat(bkg,1,1,3);
thresh=.55; % Comp
thresh=.5 % TM 

mxv=rot90(imscale(max(mxVals,thresh),256));
% rgb(:,:,1)=rgb(:,:,1)+0.2*mxv;
rgb(:,:,2)=max(rgb(:,:,2),mxv);
image(uint8(rgb))
 title(projsName,'interpreter','none');

end; % for miInd
return

%%


