% rsRefineMembraneModel3
% Newer version with Relion compatibility.
% 
% Use linear least-squares fitting to deduce a membrane profile from the
% vesicles in an image.  Given an info structure mi and a merged image, use
% the stored vesicle centers and radii to create basis functions of shells
% 1 pixel thick, which are then used in the fitting.


ds=4;           % downsampling relative to raw image
forceNewFile=0;
symmetrize=1;   % force a symmetrical profile
useAllVesicles=1;  % otherwise use the vesicles selected in SimpleRSPicker

% Modifications to the ctfs
bFactor=80;
alpha0=.02;

writeMiFile=0;

% Fitting
modelHalfwidth=50;   % angstroms
fitDC=1;        % also fit a constant offset to the image.

figure(1); clf; SetGrayscale;

% If we've already run this, we don't ask for a filename
if ~(exist('fname','var') && exist('rootPath','var')) || forceNewFile  % not already defined
    [fname pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
    % pa
    [rootPath infoPath]=ParsePath(pa);
    if ~iscell(fname)
        fname={fname};
    end;
end;

cd(rootPath);
nfiles=numel(fname);
% nfiles=1;%%%%
for fileIndex=1:numel(fname)
    disp(['Reading ' fname{fileIndex}]);
    mi=ReadMiFile([infoPath fname{fileIndex}]);
    %%
    %     iname=[mi.procPath mi.baseFilename 'm.mrc'];
    %     disp(['Reading image: ' iname]);
    m0=meLoadNormalizedImage(mi,mi.padImageSize);
    ds0=1;
    n=mi.padImageSize/ds;
    m=Downsample(m0,n);
    subplot(2,2,1);
    imacs(GaussFilt(m,256/n(1)));
    drawnow;

    if ~isfield(mi,'tempPath')
        mi.tempPath='Temp';
    end;
    if ~exist(mi.tempPath,'dir')
        mkdir(mi.tempPath);
    end

%     %     Try to read the vesicle model image if it exists
%     vname=[mi.tempPath mi.baseFilename 'v.mrc'];
%     if 0 %exist(vname,'file')
%         disp(['Reading vesicle image: ' vname]);
%         v0=ReadEMFile(vname);
%         v=Downsample(v0,n);
%     else
%         disp('Making model vesicles');
%         v=meMakeModelVesicles(mi,n);
% 
% if subtractUnfittedVesicles
%     msub=m-v;  % basic subtracted image
% end;
% drawnow;
    mres=m; % this is what we'll fit.
    msub=m;
    if useAllVesicles
        mi.vesicleModel=Crop(ones(27,1),33);
        mi.vesicle.ok= mi.vesicle.r<150 & ~isnan(mi.vesicle.s(:,1));

        numOk=sum(mi.vesicle.ok(:,3))
        subplot(2,2,2);
        imacs(GaussFilt(mres,256/n(1)));
        vesIndices=find(mi.vesicle.ok(:,3));
    else % use only the selected ones
        vesPicks=find(mi.particle.picks(:,3)==2);  % type=2 means a picked vesicle.
        vesIndices=mi.particle.picks(vesPicks,4);
        nv=numel(vesIndices);
        disp([num2str(nv) ' vesicles to be fitted']);
    end;
%     if numel(vesIndices)>0  % restore these vesicles
%         vr=meMakeModelVesicles(mi,n,vesIndices);
%         mres=msub+vr;
%         subplot(2,2,2);
%         imacs(GaussFilt(mres,256/n(1)));
%     end;
%     drawnow;
    %%

    % Create an mi structure for modification
    mi1=mi;
    %         Force all the vesicles to the same amplitude
    mi1.vesicle.s=ones(size(mi1.vesicle.s,1),1)*median(mi.vesicle.s(mi.vesicle.ok(:,1),1,1),'omitnan');
    % Make basis images
    nvs=round(modelHalfwidth/(ds*mi.pixA));
    nxs=nvs*ds;  % halfwidth of the oversampled model.
    nm=2*nvs+1;  % number of points in the model
    nx=2*nxs+1;
    if symmetrize
        nt=nvs+1;  % number of terms in the model fit
    else
        nt=nm;
    end;
    R=zeros(prod(n), nt);    % basis functions
    R0=zeros([n nt]);   % basis without CTF

    %%
    disp(['Computing ' num2str(nt) ' terms']);
    tic
    for i=1:nt % For each term, make the vesicle model vm one (or two symmetric) delta functions.
        j=(i-1-nvs)*ds+nxs+1;
        vm1=zeros(nx,1);
        vm1(j)=1;  % scale it according to thickness
        if symmetrize
            vm1(nx-j+1)=1;
        end;
        mi1.vesicleModel=vm1*ds;  % removed the mi.pixA factor
        %         mi.vesicleModel=vm1*mi.pixA;
        %         r0=meMakeModelVesicles(mi,n*ovs,vesIndices,0,0); % no filtering
        r0=meMakeModelVesicles(mi1,n,vesIndices,0,0); % no filtering
        R0(:,:,i)=r0; % ith model image

        q=Crop(r0,256);
        subplot(2,2,4);
        imacs(q);
        title(i);
        subplot(2,2,3);
        plot(sect(q));
        drawnow;
    end;
    disp('done');
    toc
    %%
    alpha0=.02;
    disp('filtering');
    for i=1:numel(mi1.ctf)
        mi1.ctf(i).alpha=alpha0;
        mi1.ctf(i).B=bFactor
        mi.ctf(i).alpha=alpha0;
    end;

    H0=ifftshift(meGetEffectiveCTF(mi1,n,ds)); % This is slow, so we compute it once.
    subplot(2,2,2);
    plot((0:n(1)/2-1)/(n(1)*ds*mi.pixA),sectr(fftshift(H0)));
    xlabel('Spatial frequency');
    ylabel('CTF')
    drawnow;
    for i=1:nt
        r=real(ifftn(fftn(R0(:,:,i)).*H0));  % Filter with the ctf
        R(:,i)=r(:);
    end;
    %     R1=reshape(R0,[prod(n) nt]);

    % Loop to fit complex CTF result
    %     disp('Fitting complex CTF');
    %     s0=mi1.vesicle.s(1);  % effective scale
    %     aScale=.5;
    niter=1;
    % Do the actual least-squares here
    mResidual=mres;
    vAccum=zeros(n);
    aAccum=zeros(nt,1);
    for j=1:niter
        if fitDC
            R(:,nt+1)=1+0*r(:);  % constant term
            a=LinLeastSquares(R,mResidual);
            a(nt+1)=[];  % remove the constant
            R(:,nt+1)=[];
        else
            a=LinLeastSquares(R,mResidual);
        end;
        %     Get the fitted model
        aAccum=aAccum+a;
        if j==1
            a0=a;
        end;
        vFitted=reshape(R*aAccum,n);
        %         v0Fitted=reshape(R1*aAccum,n);
        %         vCFitted = meSimulateComplexCTF(v0Fitted*aScale*s0,mi1)/(aScale*s0);
        %         mResidual=mres-vCFitted;
        %         mResidualC=Crop(mResidual,n/2);
        %         subplot(2,2,4);
        %         imacs(mResidualC);
        %         subplot(2,2,2);
        %         plot(sect(mResidualC));
        title(['Iteration ' num2str(j)]);
        subplot(2,2,3);
        plot([a0 aAccum]);
        %         err=mResidual(:);
        %         title(err'*err);
        drawnow;
    end;
    a=aAccum;
    %

    msub2=mres-vFitted;
    msub2c=Crop(msub2,n/2);
    msub1c=Crop(msub,n/2);
    imacs(msub2c);

    if symmetrize
        ax=[a ; flipud(a(1:nt-1))];  % full profile
    else
        ax=a;
    end;
    subplot(2,2,3);
    plot(ax,'.-','markersize',10); hold on;
    plot(0*ax,'k-'); hold off;
    a0=a;
    asym=a;
    %
    % upsample the model
    amu=meDownsampleVesicleModel(ax,1/ds);
    %

    % Try subtracting the model
    mi.vesicleModel=amu;
    %         r1=meMakeModelVesicles(mi,n,vesIndices); %  filtering on
    %
    % subplot(2,2,4);
    % msub2=mres-r1/2;
    % imacs(GaussFilt(msub2,128/n(1)));
    % drawnow;
    %
    %

    % Try to force zero baseline.
    pedpt=0;  % number of points remaining before zero crossing
    a=a0;
    pedestal=ones(nt,1);
    spike=zeros(nt,1);
    a1=a;
    if ~symmetrize
        a1(ceil((nt+1)/2):nt)=1;  % Force the 2nd half to be >0
    end;
    p1=find(a1<=0,1,'last');
    if numel(p1)<1 || p1<1+pedpt
        p1=1+pedpt;
    end;
    pedestal(1:p1-pedpt)=0;
    spike(p1-pedpt)=1;
    if symmetrize
        %         pedestal(np-p1+2:np)=0;
        %         spike(np-p1+2)=1;
    else
        a2=a;
        a2(1:ceil((nt+1)/2))=1;  % Force first half ot be positive
        p2=find(a2<=0,1,'first');
        if numel(p2)<1 || p2>nt-1-pedpt
            p2=nt-1-pedpt;
        end;
        pedestal(p2+pedpt:nt)=0;
        spike(p2+pedpt)=1;
    end;
    a=a.*pedestal;
    % R1=[R*a R*pedestal R*spike];
    R1=[R*a R*pedestal];
    a1=LinLeastSquares(R1,mres);

    % d1=[a pedestal spike]*a1;
    d1=[a pedestal]*a1;
    if symmetrize
        d1=[d1; flipud(d1(1:nt-1))];
    end;
    d1x=meDownsampleVesicleModel(d1,1/ds);
    d1x=circshift(d1x,-1);

    ndx=numel(d1x);
    ndxc=ndx-2*(ds*(p1-1)+1);
    d1xc=Crop(d1x,ndxc);
    d1xc=d1xc*2/max(d1xc);  % force maximum to be 2.
    subplot(2,2,1);
    plot(1:ndxc,d1xc,'.-',(ndxc+1)/2,d1xc((ndxc+1)/2),'r+','markersize',10);  % final composite profile
    xlabel([num2str(mi.pixA) ' A per pixel']);
    drawnow;

    outName=[mi.tempPath mi.baseFilename 'ds' num2str(ds) 'sym' num2str(symmetrize)...
        'all' num2str(useAllVesicles) 'bsc' num2str(bScale) 'binv' num2str(bInv)...
        'alph' num2str(alpha0)];
    print('-djpeg','-r300',[outName '.jpg']);  % save the plots as a jpeg.
    disp('Figure image:');
    disp(outName);

    vm.vesicleModel=d1xc;
    vm.pixA=mi.pixA;
    save([outName '.mat'],'vm');
    if writeMiFile
        mi.vesicleModel=d1xc;
        disp(['Writing ' infoPath fname{fileIndex}]);
        save([infoPath fname{fileIndex}],'mi');
    end;
end;