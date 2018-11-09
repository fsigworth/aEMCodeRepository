function meInverseFilterAuto(fname,mpars)
% Given a good subtraction of vesicles, we search for regions of the
% micrograph with low variance, and mark them as background boxes.
% Then we compute an inverse filter.
% Modified to give the absolute spectral densities in the noise model.

if nargin<1
    fname=[];
end;
if nargin<2
    mpars=struct;
end;

pars.writeMiFile=1;
pars.useUnsubImage=0; % set to 1 to get inverse filter without mbn subtraction.
pars.nb=256;
%pars.nb=64; %%%%
pars.useSmallImage=0;

pars.writeFigs=1;

pars=SetOptionValues(pars,mpars);


figNameTrunc=3;
figDir='Jpeg/';
dfc=.05;  % Gauss filtering for showing image
f0=.002; % A^1  Lower frequency cutoff
f1=.06; % A^1  Upper
%nb=64;  % box size
%nb=128;  % box size
nb=pars.nb  %%% big box size
nBoxes=14;
fHP=.005;  % assumed highpass filter
doAmpCorrection=0;
skipAbsentImages=1;

figure(1);
% Set the size so square images look square, and printing preserves the
% aspect ratio.
pos=get(gcf,'position');
set(gcf,'Position',[pos(1:2) 845 800],'PaperPositionMode','auto');
colormap(jet(256));

% Have the user select some mi files: boilerplate
%if ~exist('fname','var') || numel(fname)<1 || ~exist('doBatchProcessing','var') || ~doBatchProcessing
if nargin<1 || numel(fname)==0
    disp('Selecting mi files');
    [fname, pa]=uigetfile('*mi.txt','Select mi files','multiselect','on');
    if isnumeric(pa) % File selection cancelled
        return
    end;
    [rootPath, infoPath]=ParsePath(pa);
    if ~iscell(fname)
        fname={fname};
    end;
    cd(rootPath);
else
    infoPath='';
    rootPath='';
end;
if ~isa(fname,'cell')
    fname={fname};
end;
if ~exist(figDir,'dir')
    disp(['Making ' figDir]);
    mkdir(figDir);
end;
%%
startIndex=1;
for fileIndex=startIndex:numel(fname)
    tic;
    % if desired, create a shortened output name
    disp(fileIndex);
    name0=fname{fileIndex};
    [pa,name,ex]=fileparts(name0);
    ptr=strfind(name,'_');
    if numel(ptr)>=figNameTrunc
        outBaseName=name(1:ptr(figNameTrunc)-1);
    else
        outBaseName=name;
    end;
    % Read the mi file
    disp(['Reading ' fname{fileIndex}]);
    mi=ReadMiFile([infoPath fname{fileIndex}]);
    disp(['Defocus ' num2str(mi.ctf(1).defocus) '  Dose' ess(numel(mi.doses)) ' ' num2str(mi.doses,3)]);
    % Get a vesicle-subtracted image
    %   - first, try to get an mv.mrc or mvz.tif file
    if pars.useUnsubImage
        nameEnd='m';
    else
        nameEnd='mv';
    end;
    if pars.useSmallImage
        suffix=[nameEnd 's.mrc'];
    else
        suffix=[nameEnd,'.mrc'];
    end;
    iname=[mi.procPath mi.baseFilename suffix];
    [name,mvPresent]=CheckForImageOrZTiff(iname);
    if mvPresent
        disp(['        ' name]);
        msub=ReadEMFile(name);
        %         m=msub;
        n=size(msub);
        ds=mi.imageSize(1)/n(1);  % downsampling factor of m
        pixA=mi.pixA*ds;    % pixel size of m
        
    else % have to load the merged image and compute vesicles
        %  - read the unsubtracted file
        iname=[mi.procPath mi.baseFilename 'm.mrc'];
        [iname,ok]=CheckForImageOrZTiff(iname);
        if ok
            disp(['        ' iname]);
            m=ReadEMFile(iname);
        else
            if skipAbsentImages
                disp([iname ' not found.']);
                continue;
            end;
            [iname,pa]=uigetfile('*m*','Find the merged image');
            cd(pa);
            m=ReadEMFile(iname);
        end;
        
        %%
        
        % Get image and pixel sizes
        n=size(m);
        ds=mi.imageSize(1)/n(1);  % downsampling factor of m
        pixA=mi.pixA*ds;    % pixel size of m
        
        mysubplot(2,2,1); % first display without subtraction
        ShowImageAndBoxes(GaussFilt(m,dfc));
        title(['doses: ' num2str(mi.doses) '  weights: ' num2str(mi.weights)]);
        drawnow;
        
        if isfield(mi,'vesicle') && numel(mi.vesicle.x)>0
            
            % get the vesicle subtraction
            %         disp('Subtracting vesicles');
            vfile=[mi.tempPath mi.baseFilename 'v.mrc'];
            ok=0;
            %     if isfield(mi,'tempPath') && exist(vfile,'file')
            %         vs=ReadMRC(vfile);
            %         ok=all(size(vs)==n);
            %     end;
            %     if ~ok
            vs=meMakeModelVesicles(mi,n);
            %     end;
            if doAmpCorrection
                %   Do a least-squares amp correction
                vesicleAmpCorrection=(vs(:)'*m(:))/(vs(:)'*vs(:))
            else
                vesicleAmpCorrection=1;
            end;
            vs=vs*vesicleAmpCorrection;
            msub=m-vs;
            mvName=[mi.procPath mi.baseFilename 'mvz.tif'];
            WriteZTiff(msub,pixA,mvName);
            disp(['written: ' mvName]);
        else
            disp('No vesicle data available, using unsub image.');
            msub=m;
        end;
    end; % if mvPresent
    
    %         Ti=meGetNoiseWhiteningFilter(mi,n);
    %     msub=real(ifftn(fftn(msub).*ifftshift(Ti)));
    
    
    mysubplot(2,2,1); % Show subtracted image
    ShowImageAndBoxes(GaussFilt(msub,dfc));
    title(['doses: ' num2str(mi.doses) '  weights: ' num2str(mi.weights) 'ds: ' num2str(ds)]);
    drawnow;
    
%     % Get the effective CTF from the merging.
%     H=ifftshift(meGetEffectiveCTF(mi,n,ds));
    
    %% Put boxes in the lowest-variance regions and compute 1D spectrum
    % -------Get the variance map--------
    %     disp('Computing the variance map');
    %     1.  Filter the image
    df=1./(n*pixA);  % frequency step in the whole image
    % create a bandpass filter from f0 to f1
    filt=fuzzymask(n,2,f1./df,f1/(4*df(1)))-fuzzymask(n,2,f0./df,f0/(4*df(1)));
    mfilt=real(ifftn(fftn(msub).*ifftshift(filt)));
    %     2.  Get the convolution box
    box=Crop(SquareWindow(nb+2,2),n);
    boxo=ifftshift(box);
    localVar=real(ifftn(fftn(mfilt.^2).*fftn(boxo)));
    mask=meGetMask(mi,n);
    mask=Crop(Crop(mask,n-ceil(max(64/ds+nb/2,n/6))),n);  % Force a band of zeros at edges
    mxVar=max(localVar(:));
    lVar=localVar.*mask+mxVar*(1-mask);
    %     subplot(2,2,3);
    %     imacs(lVar)
    drawnow;
    %     3. Get some boxes
    lv2=lVar;
    blankBox=ifftshift(Crop(ones(nb*2),n));
    boxX=zeros(1,nBoxes);
    boxY=zeros(1,nBoxes);
    imgs=zeros(nb,nb,nBoxes);
    for i=1:nBoxes
        [val ix iy]=max2d(-lv2);
        lv2=max(lv2,circshift(blankBox,[ix iy])*mxVar);  % blank the region
        imgs(:,:,i)=ExtractImage(msub,[ix iy],nb);
        %         disp([-val ix iy]);
        boxX(i)=(ix-1)*ds;
        boxY(i)=(iy-1)*ds;
    end;
    mysubplot(2,2,1);
    ShowImageAndBoxes(GaussFilt(msub,dfc),[boxX' boxY']/ds,nb,8,[1 1 0]);
    title(['doses: ' num2str(mi.doses) '  weights: ' num2str(mi.weights) '  ds: ' num2str(ds)]);
    drawnow;
    
    %        4.  Compute average spectrum
    sp=mean(RadialPowerSpectrum(imgs,1),2)*pixA^2;  % To get the general spectral
    %     density in units of A^2 because we multiply this by pixA^2
    dfb=1/(pixA*nb); % frequency step in boxes
    freqs=(0:nb/2-1)'*dfb; % frequencies in spectrum
    c=sectr(meGetEffectiveCTF(mi,nb,ds));  % corresponding ctf
    
    subplot(4,2,6); % Show the spectrum and the ctf
    plot(freqs,[sp/40 c]);
    
    %% Fit the power spectrum with NoiseModel.*c^2+shot
    noiseModelFcn='NoiseModel1';
    niters=1000;
    
    % Get the effective CTFs
    effctf=sectr(meGetEffectiveCTF(mi,nb,ds));
    %     effctfx=sectr(meGetEffectiveCTF(mi,nb*10,ds));
    
    % fit the noise model with multiple families of parameters
    %     nPSets=1;
    nPSets=3;  % try this many sets of parameters
    % Initialize parameters:
    q=sp(nb/4);  % get a rough sample of the spectral density
    
    % noise model 1 is
    %     spec=(ag*gauss(sigma)+af1).*(.01./f).^f1exp+af2*(.01./f).^f2exp;
    %     shot=s0./(1+(f/bf).^2);
    %     Really simple model, with a/f^n+b for spect, c for shot.
    %        af1 af2  ag   sigma  bf  s0 f1exp f2exp
    p0=[     3*q  0    0     1   100   q  2    1 ];
    ac=[      1   0    0     0    0    1  0    0 ];
    
    %        af1 af2 ag    sigma  bf  s0 f1exp f2exp
    p0(2,:)=[ 1*q 1*q 10*q .007  100  q   0.5   1.5];
    ac(2,:)=[ 1   1    1    1    0   1    1     1 ];
    %        af1 af2 ag    sigma  bf  s0 f1exp f2exp
    %     p0(2,:)=[ 1*q 1*q 10*q .007  100  q   1.5   1.5];
    %     ac(2,:)=[ 1   1    1    1    0   1    1     1 ];
    %
    p0(3,:)=[1*q 1*q  0*q .007  100  q   1.5   1.5];
    ac(3,:)=[1   1    0    0    0    1    1     1 ];
    
    ps=p0(1:nPSets,:);   % Store the final parameters
    errs=zeros(nPSets,1); % store the final errors.
    
    for jp=1:nPSets;  % loop over sets of parameters
        p=Simplex('init',p0(jp,:),ac(jp,:));
        for i=1:niters
            [spec, shot]=eval([noiseModelFcn '(freqs,p)']);
            model=spec.*effctf.^2+shot;
            d=(model-sp);
            err=d'*d;
            err=err+1e6*sum(p<0);  % penalty for neg. values
            p=Simplex(err,p);
            
            if mod(i,100)==0
                subplot(4,2,6);
                plot(freqs,sp,'k.',freqs,model,'b-');
                title(['model: ' num2str(jp) '  iter: ' num2str(i)]);
                xlabel('Spatial frequency, A^{-1}');
                ylabel('Spectral density, A^2')
                subplot(2,2,3);
                semilogy(freqs,[spec shot]);
                ylabel('Model: spec, shot density, A^2');
                xlabel('Spatial frequency, A^{-1}');
                drawnow;
            end;
        end;
        errs(jp)=err;
        ps(jp,:)=Simplex('centroid');
    end;
    errs
    %     Pick the one which converged better
    [~, jBest]=min(errs);
    p=ps(jBest,:);
    
    disp('      af1       af2       ag      sigma       bf        s0      f1exp     f2exp');
    disp(p);
    
    % Show the final results
    [spec, shot]=eval([noiseModelFcn '(freqs,p)']);
    model=spec.*effctf.^2+shot;
    
    subplot(4,2,6);
    plot(freqs,sp,'k.',freqs,model,'b-');
    title(['model: ' num2str(jBest)]);
    xlabel('Spatial frequency, A^{-1}');
    ylabel('Spectral density, �^2');
    
    
    subplot(2,2,3);
    semilogy(freqs,[spec shot model sp]);
    ylabel('Model: spec, shot density');
    xlabel('Spatial frequency, A^{-1}');
    title(['Shot density: ' num2str(mean(shot)) ' �^2']);
    drawnow;
    
    mi=meStoreNoiseModel(p,noiseModelFcn,mi);
    
    if pars.writeMiFile
        mi.log{end+1,1}=['meInverseFilterAuto ' TimeStamp];
        WriteMiFile(mi,[rootPath infoPath fname{fileIndex}]);
        disp(['Updated ' fname{fileIndex}]);
    else
        disp('mi file not updated.');
    end;
    
    % Whiten a binned image for display
    dsDisplay=2;
    n2=n/dsDisplay;
    ds2=ds*dsDisplay;
    msub2=BinImage(msub,dsDisplay);
    pixA2=pixA*dsDisplay;
    
    Ti=meGetNoiseWhiteningFilter(mi,n2,ds2,1,fHP); % (1 means 1 zero)
    f2d=RadiusNorm(n2)/pixA2;
        subplot(4,2,8);                         %%
    plot(sectr(f2d),sectr(Ti));
    axis([0 inf 0 1.2]);
    %     axis([0 inf 0 1.02*max(sectr(Ti))]);    %%
    title('Prewhitening filter');           %%
    xlabel('Spatial frequency');    %%
    ylabel('T(s)');
    drawnow;
    msubf=real(ifftn(fftn(msub2).*ifftshift(Ti)));
    %%
    mysubplot(2,2,2);
    %     ShowImageAndBoxes(GaussFilt(msubf,.1),[boxX' boxY']/ds,nb,8,[1 1 0]);
    ShowImageAndBoxes(GaussFilt(msubf,dfc*2));
    axis off;
    %     disp(mi.baseFilename);
    title(['Filtered: ' mi.baseFilename],'interpreter','none');
    drawnow;
    %     subplot(2,2,3);
    %     plot(sectr(f2d)/pixA,[RadialPowerSpectrum(msubf) RadialPowerSpectrum(msub)]);
    %     title(num2str([p(1) p(6)]));
    %     drawnow;
    
    if pars.writeFigs
        figName=[figDir outBaseName 'pwf.jpg'];
        print('-djpeg','-r150',figName);  % save the CTF window.
        disp(['Wrote figure: ' figName]);
    end;
    
    %     if writeOut
    %         outName=[mi.baseFilename 'mw.mrc'];
    %         disp(['Writing ' outName]);
    %         WriteMRC(mf,pixA,[mi.procPath outName]);
    %     end;
    toc
    disp(' ');
end;