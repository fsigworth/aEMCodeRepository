% spPickingPreprocessor4
% This program does the correlations for SimpleRSPicker on conventional
% single-particle data.  It writes a file
% *rscc.mat which contains the information.  It also places
% Some of the 2D arrays written out in the *rscc.mat file.
% mxCC is the computed amplitude of the putative particle at each pixel
% position.  The amplitude is relative to the
%  amplitude of the reference.
%  amxTemplInds gives the index of the reference giving that value. This is the
%  index k into eigenSet.vList(:,k).
% mxVesInds gives the best-match or else the nearest vesicle to each pixel.
% mxRsos is 1 for a right-side-out particle at that position.  It
% is zero for inside-out and is also zero outside vesicle regions.
% mxVars is the local variance, computed as thesquare of the cc of
% the image with the first 1 or a few terms of the eigenimage expansion.

kvMode=2;  % otherwise, AMPAR
useNoiseWhitening=1;
pwFiltPars=[.002 0; .01 .5];  % generic pw filter parameters
pwFiltPars=[.002 0];  % generic pw filter parameters
fHP=.001; % highpass in A^-1
overrideDsm=2;  % downsample the merged image by this factor.
overrideB=40;
jetScale=600;
localVarRadius=200;
membraneThicknessA=60;
maskPaddingA=200;
fakeImageIndex=0;
nterms=30;
deleteNWF=1;

s0=1e-3;  % amplitude contrast factor
ppVals=struct;
switch kvMode
    case 2
        ppVals.membraneOffsetA=0;
        mapName='TrpBin2.mrc';
        dimShift=0;
        symmetry=4;
    case 1
        ppVals.membraneOffsetA = 52;  % membrane is 52 Å above particle center for RSO particle
        %    Hence a particle center > r-52Å implies an inside-out particle.
        %     In turn, the center of beta subunits is about 50 Å below the particle
        %     center
        %
        %     mapName='/Users/fred/Structures/kv1.2/KvMap.mrc';
        %     mapName='~/Structures/kv1.2/KvMap.mrc'; % weak membrane subtraction
        mapName='~/Structures/kv1.2/KvMapMbnSub.mrc';  % stronger membrane subtraction.
        dimShift=2;
        symmetry=4;
    otherwise  % AMPAR
        ppVals.membraneOffsetA = -70;  % Membrane center is this distance from particle center
        mapName='/Users/fred/Structures/AMPAR/3KG2map58.mrc';
        dimShift=0;
        symmetry=2;
end;
% ppVals.nAlpha=32; % about 10 degrees
% ppVals.nBeta=12;  % even is best; number of betas per hemisphere.
% ppVals.nGamma=16/symmetry;

% new tighter values
ppVals.nAlpha=18;
ppVals.nBeta=6;  % even is best; number of betas per hemisphere.
ppVals.nGamma=16/symmetry;
% This struct is copied to the mi file in rspLoadFiles.

% outputImageSize=1024;  % size of output images.
outputImageSize=960;  % size of output images.
% outputImageSize=768;
% nterms=1;
simulateImage=0;

doRelVarFigure=0;  % special figure of variance vs. number of terms
% nterms=40;       % use more terms for this figure.


% mapName='/Volumes/TetraData/Structures/AMPAR/3KG2mapsub5.8A.mrc';

% Have the user select some mi files: boilerplate
if ~exist('fname','var') || ~exist('doBatchProcessing','var') || ~doBatchProcessing
    [fname, pa]=uigetfile('*mi.txt','Select mi files','multiselect','on');
    if isnumeric(pa) % File selection cancelled
        return
    end;
    [rootPath, infoPath]=ParsePath(pa);
    if ~iscell(fname)
        fname={fname};
    end;
    cd(rootPath);
end;

oldPixA=0;
%%
for fileIndex=1:numel(fname); % Operate on a single micrograph
    tic
    miName=[infoPath fname{fileIndex}];
    disp(['Reading ' miName]);
    mi=ReadMiFile(miName);
    mi.basePath=rootPath;
    %     Pick up the merged image
    [origImg, mergeFullPath]=meReadMergedImage(mi);
    %%
    %     Possibly downsample it
    n0=size(origImg);
    dsm=ceil(n0(1)/outputImageSize); % Downsampling factor from merged image
    if overrideDsm
        dsm=overrideDsm;
    end;
    n=n0(1)/dsm;
    if dsm~=1  % further downsampling
        m1=Downsample(origImg,n);
    else
        m1=origImg;
    end;
    n=size(m1);
    ctr=n/2+1;  % n must be even
    ds=mi.imageSize(1)/n(1);  % downsampling relative to original micrograph
    pixA=mi.pixA*ds;
    
    %         useNoiseWhitening=(numel(mi.noiseModelPars)>0);
    if  useNoiseWhitening
        disp('Noise whitening');
        [H,ok]=meGetNoiseWhiteningFilter(mi,n,0,1,fHP*pixA);
        if ~ok || deleteNWF
            mi.noiseModelPars=[];
            if ok && deleteNWF
                disp('Deleting the noise-whitening filter pars');
                WriteMiFile(mi,miName);
            end;
            disp('Using a generic filter');
            freqs=RadiusNorm(n)/pixA+1e-6;  % in inverse Å, prevent divide by zero
            H=ones(n,'single');
            for i=1:size(pwFiltPars,1)  % product of (gauss + const) terms.
                f0=pwFiltPars(i,1);
                a=pwFiltPars(i,2);
                h=exp(-(f0./freqs).^2);
                H=H.*(a+(1-a)*h);
            end;
        end;
        m2=real(ifftn(fftn(m1).*ifftshift(H)));
    else
        disp('No specimen-noise whitening');
        m2=m1;
        H=1;
    end;
    %
    figure(1);
    colormap jet(256);
    mysubplot(221);
    imaga(imscale(GaussFilt(m2,.3),256,.001));
    axis off;
    title(['Prewhitened: ' fname{fileIndex}],'interpreter','none');
    drawnow;
    %%
    if pixA~=oldPixA  % We haven't already made templates of the correct size
        oldPixA=pixA;
        % Load the 3D map
        disp('Loading the 3D map');
        [origMap, mpixA]=ReadEMFile(mapName);
        % mpixA=1.8 %%%%%%%%
        nt1=size(origMap,1)*mpixA/pixA;  % final effective map size
        nt=ceil(nt1/8)*8;
        [map, finalmag]=DownsampleGeneral(origMap,nt,mpixA/pixA);
        map=map*pixA*s0;  % approx amplitude correction (V-A scaling) x sigma
        map=shiftdim(map,dimShift);
        % magnifications=[mpixA/pixA finalmag]
        
        %% Create the list of angles for the templates
        membraneOffset=ppVals.membraneOffsetA/pixA;
        membraneOffsetA=ppVals.membraneOffsetA;
        gammaStep=360/(symmetry*ppVals.nGamma);
        
        %         hemiAngles run from alpha=[0..360) and beta=[0..90)
        [hemiAngles, angleInds]=rsListHemisphereAngles(ppVals.nAlpha, ppVals.nBeta);
        nHemiAngles=size(hemiAngles,1);  % list of [alpha beta]
        nHemi=2;  % both hemispheres
        
        angleList=zeros(nHemi,ppVals.nGamma,nHemiAngles,3,'single');
        for j=1:ppVals.nGamma;
            gamma=(j-1)*gammaStep;
            for k=1:nHemiAngles
                angleList(1,j,k,:)=[hemiAngles(k,:) gamma];
                angleList(2,j,k,:)=[hemiAngles(k,1) 180-hemiAngles(k,2) gamma];
                %                 angleList(2,j,k,:)=[[0 180]-hemiAngles(k,:) gamma];
            end;
        end;
        nAngles=numel(angleList)/3;
        
        % angle list is of size
        % (nHemi x nGamma x nHemiAngles, 3)  where nHemi=2 is the
        % number of hemispheres, and the hemisphere index is
        % fastest-varying.  The rows are [alpha beta gamma]. Thus
        % the beta angles alternate such that, if betastep=1, they are
        % are (0, 180) for each gamma, then (1, 179) for each gamma, up to 89, 91.
        % that is, the same projected position is described twice.
        
        %% Make the templates
        
        disp(['Making ' num2str(nAngles) ' templates']);
        
        % allTemplates=rsMakeTemplatesQuick(angleList,map);
        allTemplates=rsMakeTemplates(reshape(angleList,nAngles,3),map);
        %             toc
        allTemplates=reshape(allTemplates,nt,nt,nHemi,ppVals.nGamma,nHemiAngles);
        
    end;  % if pixA
    
    %% Filter the templates according to the CTF
    [nt, nt, nHemi, ppVals.nGamma, nHemiAngles]=size(allTemplates);
    nAngles=nHemi*ppVals.nGamma*nHemiAngles;
    ne=NextNiceNumber(nt*1.3);  % increase the size to allow CTF rings
    mi1=mi;
    if overrideB
        for i=1:numel(mi.ctf)
            mi.ctf(i).B=overrideB;
        end;
    end;
    ctf=meGetEffectiveCTF(mi1,ne,ds);  % put in dqe, pw filter.
    if useNoiseWhitening % make a pre-whitening filter for the references
        pwfRef=meGetNoiseWhiteningFilter(mi1,ne,ds);
        if numel(pwfRef)<2
            disp('Using a generic filter');
            freqs=RadiusNorm(ne)/pixA+1e-6;  % in inverse Å, prevent divide by zero
            pwfRef=ones(ne,ne,'single');
            for i=1:size(pwFiltPars,1)  % product of (gauss + const) terms.
                f0=pwFiltPars(i,1);
                a=pwFiltPars(i,2);
                h=exp(-(f0./freqs).^2);
                pwfRef=pwfRef.*(a+(1-a)*h);
            end;
        end;
    else
        pwfRef=ones(ne,ne,'single');
    end;
    ctf=ctf.*pwfRef;
    
    % Pad the templates to avoid ctf artifacts
    xTemplates=Crop(reshape(allTemplates,nt,nt,nAngles),ne,1);
    nim=size(xTemplates,3);
    %     operate with the CTF and mask
    H=ifftshift(ctf);
    msk=fuzzymask(ne,2,0.45*ne,.1*ne);
    for i=1:nim
        xTemplates(:,:,i)=real(ifftn(fftn(xTemplates(:,:,i)).*H)).*msk;
    end;
    xTemplates=reshape(xTemplates,ne,ne,nHemi,ppVals.nGamma,nHemiAngles);
    ntstr=num2str(nt);
    nestr=num2str(ne);
    disp(['Templates expanded from ' ntstr 'x' ntstr ' to ' nestr 'x' nestr]);
    
    %% Make the eigenreferences
    eigenSet=rsMakeEigenrefs(xTemplates,nterms);
    %%
    mic=mi;
    
    
    %     makeFakeImage=ceil(rand*nRefs);
    
    if fakeImageIndex
        sigmaN=.15;
        fakeImg=-Crop(xTemplates(:,:,fakeImageIndex),n); % scale up one reference
        m2=sigmaN*randn(n)+fakeImg;
        m1=Downsample(m2,n0);
        mysubplot(221);
        imags(Crop(fakeImg,128));
        mysubplot(223);
        imags(Crop(m2,128));
    end;
    
    figure(1);
    mysubplot(222);
    %     imaga(imscale(GaussFilt(m1,.05),256,.0005));
    imags(GaussFilt(m2,.1));
    axis off;
    drawnow;
    
    globalMask=meGetMask(mi,n);
    
    %
    
    %  ----------Evaluate the cc ----------
    disp('Evaluating cross-correlations');
    disp([' using ' num2str(nterms) ' terms']);
    
    mxVesInds=uint16(zeros(n));
    mxVars=ones(n,'single');
    mxRsos=ones(n,'single');
    mxDist=single(ones(n))*max(n);
    nulls=zeros(n,'single');
    badVesMask=false(n);
    msklRadius=round(localVarRadius/pixA);
    membraneThicknessPix=ceil(membraneThicknessA/pixA);
    nl=2*ceil(1.2*msklRadius);
    mskl=fuzzymask(nl,2,msklRadius,0.12*msklRadius);  % local mask
    npts=sum(mskl(:));
    
    maskPadding=ceil(maskPaddingA/pixA);
    
    %                 % ---------- Get the correlation function here ------------
    [gMaxVals, mxTemplInds, mxNCC]=spPickingCorrelation6(-m2,mic,eigenSet);
    
    mysubplot(223);
    imacs(gMaxVals);  % Show the CC function
    axis off;
    drawnow;
    
    mysubplot(224);
    plot(sect(gMaxVals));
    title('NCC');
    
    mxCC=gMaxVals;
    eigenImgs=eigenSet.imgs;
    vList=single(eigenSet.vList);
    mVesGood=zeros(n,'single');
    mVesBad=zeros(n,'single');
    %         Quantize the vesicle models
    % mVesGood=.002*round(500*mVesGood);
    % mVesBad=.002*round(500*mVesBad);
    log=['spPickingPreprocessor ' TimeStamp];
    spMode=1;
    
    partRadius=18; % ????
    save([mi.procPath mi.baseFilename 'rscc.mat'],'mxCC','mxVars','mxVesInds',...
        'mxDist','mxTemplInds','mxRsos','partRadius', 'membraneOffsetA','ds',...
        'badVesMask','eigenImgs','vList','angleList','ppVals','pwfRef',...
        'mVesGood','mVesBad','spMode','log');
    disp(['written: ' mi.procPath mi.baseFilename 'rscc.mat']);
    toc
    disp(' ');
    
    %
    if fakeImageIndex
        figure(2);
        ip=Crop(mxTemplInds,5)
        ic=round(100*Crop(mxCC,5))
        inc=round(100*Crop(mxNCC,5))
        [val,ptr]=max(ic(:));
        iTempl=ip(ptr)
        
        subplot(231);
        imags(xTemplates(:,:,fakeImageIndex));
        
        subplot(234);
        imags(xTemplates(:,:,iTempl));
        
        eis=reshape(eigenSet.imgs,ne^2,nterms);
        evs=reshape(eigenSet.vList,nterms,nRefs);
        
        recon=reshape(eis*evs(:,fakeImageIndex),ne,ne);
        subplot(232);
        imags(recon);
        title(['true: ' num2str(fakeImageIndex)]);
        
        recon=reshape(eis*evs(:,iTempl),ne,ne);
        subplot(235);
        imags(recon);
        title(iTempl);
        
        subplot(233);
        imags(Crop(mxNCC,32));
        hold on;
        plot(17,17,'r+');
        hold off;
        title('NCC');
        
        subplot(236);
        imags(Crop(mxCC,32));
        hold on;
        plot(17,17,'r+');
        hold off;
        title(['Amp= ' num2str(val)]);
        figure(1);
        imags(GaussFilt(Crop(m2,128),.2));
        
    end;
end;

%%
