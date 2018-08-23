% rsPickingPreprocessor

membraneOffsetA = -70;  % Membrane center is this distance from particle center
localVarRadius=100;  % angstroms
maskPaddingA=5;

ds1=2;  % extra downsampling compared to the merged image.
showTemplates=0;
simulateImage=0;

mapName='/Volumes/TetraData/Structures/AMPAR/3KG2map58.mrc';
% mapName='/Volumes/TetraData/Structures/AMPAR/3KG2mapsub5.8A.mrc';

% Have the user select some mi files
[fname pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
[rootPath infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;
cd(rootPath);

oldPixA=0;

for fileIndex=1:numel(fname); % Operate on a single micrograph
    miName=[infoPath fname{fileIndex}];
    disp(['Reading ' miName]);
    load(miName);
    if ~isfield(mi,'vesicle.ok')
        amps=mi.vesicle.s;
        minS=0;
        maxS=inf;
        mi.vesicle.ok=(amps>minS) & (amps<maxS);
    end;
    
%     Pick up the merged image
    imgExts={'mc.mrc','m.mrc'};
    i=0;
    ok=0;
    while i<=numel(imgExts)&&~ok
        i=i+1;
        imgName=[mi.procPath mi.baseFilename imgExts{i}];
        if exist(imgName,'file')
            ok=1;
            origImg=ReadEMFile(imgName);
        end;
    end;
    if ~ok
        error(['Couldn''t open the image ' imgName]);
    end;

    %     Possibly downsample it
    n=size(origImg);    
    if ds1>1  % further downsampling
        n=n/ds1;
        m0=Downsample(origImg,n);
    else
        m0=origImg;
    end;
    
    ds=mi.imageSize(1)/n(1);
    pixA=mi.pixA*ds;
    
    disp('Making model vesicles');
    ves=meMakeModelVesicles(mi,n);
    m1=m0-ves;
    
    useNoiseWhitening=(numel(mi.noiseModelPars)>0);
    if  useNoiseWhitening
        disp('Noise whitening');
        H=meGetNoiseWhiteningFilter(mi,size(m1));
        m=real(ifftn(fftn(m1).*ifftshift(H)));
    else
        disp('No specimen-noise whitening');
        m=m1;
    end;
    
    figure(1); clf;
    SetGrayscale;
    imac(imscale(GaussFilt(m,.2),256,.001));
    drawnow;

    if pixA~=oldPixA  % We haven't already made templates of the correct size
        % Load the 3D map
        disp('Loading the 3D map');
        [origMap s]=ReadMRC(mapName);
        mpixA=s.pixA;
        nt1=size(origMap,1)*mpixA/pixA;  % final effective map size
        nt=ceil(nt1/8)*8;
        [map finalmag]=DownsampleGeneral(origMap,nt,mpixA/pixA);
        map=map*pixA;  % approx amplitude correction (V-A scaling)

        % magnifications=[mpixA/pixA finalmag]

        %% Create the list of angles for the templates
        ds=2*ds1;  % We started with a 2x downsampled merged image.  Overall downsampling is 4.
        membraneOffset=membraneOffsetA/pixA;
        nAlpha=32; % about 10 degrees
        nBeta=12;  % even is best.  Number on hemisphere.
        nGamma=8;
        symmetry=2;
        gammaStep=360/(symmetry*nGamma);
        nterms=28;
    %     nterms=12;

        [hemiAngles angleInds]=rsListHemisphereAngles(nAlpha, nBeta);
        nHemiAngles=size(hemiAngles,1);
        nHemi=2;  % both hemispheres

        angleList=zeros(nHemi,nGamma,nHemiAngles,3);
        for j=1:nGamma;
            gamma=(j-1)*gammaStep;
            for k=1:nHemiAngles
                angleList(1,j,k,:)=[hemiAngles(k,:) gamma];
                angleList(2,j,k,:)=[[0 180]-hemiAngles(k,:) gamma];
            end;
        end;
        nAngles=numel(angleList)/3;

        % angle list is of size
        % (nHemi x nGamma x nHemiAngles, 3)  where nHemi=2 is the number of hemispheres.
        % Note that the beta angles alternate such that, if betastep=1, they are
        % are (0, 180) for each gamma, then (1, 179) for each gamma, up to 89, 91.
        % that is, the same projected position is described twice.

        %% Make the templates

        disp(['Making ' num2str(nAngles) ' templates']);

        tic
        % allTemplates=rsMakeTemplatesQuick(angleList,map);
        allTemplates=rsMakeTemplates(reshape(angleList,nAngles,3),map);
        toc

        % angleList=reshape(angleList,nHemiAngles,2,nGamma,3);
        allTemplates=reshape(allTemplates,nt,nt,nHemi,nGamma,nHemiAngles);
    
    end;
    
    %% Filter the templates according to the CTF
    [nt nt nHemi nGamma nHemiAngles]=size(allTemplates);
    nAngles=nHemi*nGamma*nHemiAngles;
    ne=NextNiceNumber(nt*1.2);  % increase the size to allow CTF rings
    ctf=meGetEffectiveCTF(mi,ne,ds);  % put in dqe, pw filter.
    if useNoiseWhitening
        ctf=ctf.*meGetNoiseWhiteningFilter(mi,ne,ds);
    end;
%     % evaluate a generic inverse filter
%     f=RadiusNorm(ne)/(mi.pixA*ds);  % frequencies for evaluating the CTF
%     hinv=a1./(1+(f1./f).^2)+(1-a1)*1./(1+(f0./f).^ex);
%     hinv(f==0)=0;
    % Pad the templates to avoid ctf artifacts
    % eigenSet=rsMakeEigenrefs(allTemplates,nterms,ctf.*hinv);
    xTemplates=Crop(reshape(allTemplates,nt,nt,nAngles),ne,1);
    nim=size(xTemplates,3);
    H=ifftshift(ctf);
    msk=fuzzymask(ne,2,0.45*ne,.1*ne);
    for i=1:nim
        xTemplates(:,:,i)=real(ifftn(fftn(xTemplates(:,:,i)).*H)).*msk;
    end;
    xTemplates=reshape(xTemplates,ne,ne,nHemi,nGamma,nHemiAngles);
    ntstr=num2str(nt);
    nestr=num2str(ne);
    disp(['Templates expanded from ' ntstr 'x' ntstr ' to ' nestr 'x' nestr]);
    
    %% Make the eigenreferences
    eigenSet=rsMakeEigenrefs(xTemplates,nterms);
    figure(2);
    SetGrayscale;
    plot(1-eigenSet.termVar);
    % figure(1);
    % ImagicDisplay(eigenSet.imgs,2);
    
    
    if showTemplates
        %%  % Show the templates and the reconstructions
        timgs=reshape(eigenSet.imgs,ne*ne,nterms);
        nG=1;
        nH=2;
        nAngs=nHemiAngles;
        rImg=single(zeros(ne,ne,2*nH,nG,nAngs));
        for k=1:nHemiAngles
            for j=1:nG
                for i=1:nH
                    %             rImg(:,:,2*i-1,j,k)=reshape(timgs*eigenSet.vList(:,i,j,k),ne,ne).*squeeze(eigenSet.ampList(i,j,k));
                    rImg(:,:,2*i-1,j,k)=reshape(timgs*eigenSet.vList(:,i,j,k),ne,ne);
                    rImg(:,:,2*i,j,k)=xTemplates(:,:,i,j,k);
                end;
            end;
        end;
        figure(1);
        ImagicDisplay2(rImg,2);
        
        %%  % Make a complete set of reconstructions for comparison
        %   and make the average power spectrum
        nTotal=nAngs*nGamma*nHemi;
        timgs=reshape(eigenSet.imgs,ne*ne,nterms);
        rImg=single(zeros(ne,ne,nTotal));
        vL=reshape(eigenSet.vList,nterms,nTotal);
        sp=zeros(ne,ne);
        sp0=zeros(ne,ne);
        for i=1:nTotal
            %     img=reshape(timgs*vL(:,i),ne,ne)*eigenSet.ampList(i);
            img=reshape(timgs*vL(:,i),ne,ne);
            sp0=sp0+abs(fftn(xTemplates(:,:,i))).^2;
            sp=sp+abs(fftn(img)).^2;
            rImg(:,:,i)=img;
        end;
        spr=Radial(fftshift(sp/nTotal));
        spr0=Radial(fftshift(sp0/nTotal));
        figure(2);
        semilogy([spr0 spr]);
        figure(3);
        plot([cumsum(spr0) cumsum(spr)]);
        ImagicDisplay2(rImg);
    end;
    % mic=rsSortVesicles(mi);  % make a copy with the vesicles sorted by position.
    mic=mi;
    
    %%  Evaluate the cc for each vesicle
    disp('Evaluating cross-correlations');
    figure(1); clf; SetGrayscale;
    imac(imscale(GaussHP(GaussFilt(origImg,.1),.005),256,.003));
    namp=1.5;
    if simulateImage
        m2=m+namp*randn(n);
    else
        m2=m;
    end;
    
    
    
    % mic.vesicle.ok=(mic.vesicle.s>minS) & (mic.vesicle.s<maxS);
    % % mic.vesicle.s=mic.vesicle.s*0+.05;  % make the amplitudes equal.
    vesInds=find(mic.vesicle.ok);
    
    disp(['Total vesicles: ' num2str(numel(vesInds))]);
    mxCC=single(zeros(n));
    mxTemplInds=uint16(zeros(n));
    mxVesInds=uint16(zeros(n));
    mxVars=single(zeros(n));
    mxRsos=single(zeros(n));
    figure(3);
    SetGrayscale;
    
    msklRadius=round(localVarRadius/pixA);
    nl=2*ceil(1.2*msklRadius);
    mskl=fuzzymask(nl,2,msklRadius,0.12*msklRadius);  % local mask
    npts=sum(mskl(:));
    
    nves=numel(vesInds);
    maskPadding=ceil(maskPaddingA/pixA);
    for j=1:nves
        i=vesInds(j);
        maskRadii=mic.vesicle.r(i)/ds+membraneOffset*[-1 1 0]+maskPadding;
        maskRadii(3)=maskRadii(3)+abs(membraneOffset)+msklRadius+maskPadding;
        % Get the single-vesicle correlation function
        [mxVals mxInds mxNCC mxValsV mxRso]=rsVesicleCorrelation4(-m2,mic,i,...
            membraneOffset,maskRadii,angleInds,eigenSet);
        ctr=round([mic.vesicle.x(i) mic.vesicle.y(i)]/ds+1);  % shift to 1-based coords
        nv=size(mxVals,1);
        
        % Do a local averaging of the squared CC in the whole vesicle.
        var=mxVals.^2;
        h=ifftshift(Crop(mskl,nv));  % average over the local mask
        filtVar=real(ifftn(fftn(var).*fftn(h))).*fuzzymask(nv,2,max(maskRadii(1:2)),.5);
        
        %     Pad to the full-sized image
        xVals=ExtractImage(mxVals,ctr,n,1);
        xTemplInds=ExtractImage(mxInds,ctr,n,1);
        xVars=ExtractImage(filtVar,ctr,n,1);
        xNCC=ExtractImage(mxNCC,ctr,n,1);
        xRso=ExtractImage(mxRso, ctr,n,1);
        
        %       Incorporate into the composite image
         q=xVals>mxCC;  % find maximum cc values
 %       q=xNCC>mxCC;  % find maximum cc values
         mxCC(q)=xVals(q); % update the values where these ones are higher.
 %       mxCC(q)=xNCC(q); % update the values where these ones are higher.
        mxVesInds(q)=i;  % the related vesicle indices
        mxTemplInds(q)=xTemplInds(q);  % the template indices
        mxRsos(q)=xRso(q); % the right-side-out flags.
        q1=xVars>mxVars;  % find maximum variance values
        mxVars(q1)=xVars(q1);
        
        if mod(i,10)==0
            figure(3);
            imacs(mxCC);  % Show the CC function
            %         imacs(mxVars);  % Show the local variances
            title(i);
            drawnow;
        end;
   
    end;  % loop over vesicles
    
    partRadius=18;
    save([mi.procPath mi.baseFilename 'rscc.mat'],'mxCC','mxVars','mxVesInds',...
        'mxTemplInds','mxRsos','m0','m1','partRadius','ds');
    disp(['written: ' mi.procPath mi.baseFilename 'rscc.mat']);
    disp(' ');
end;
