% rsRefineMembraneModel3 -- experimental version
% Fit the differently-defocused images separately.
%
% Use linear least-squares fitting to deduce a membrane profile from the
% vesicles in an image.  Given an info structure mi and a merged image, use
% the stored vesicle centers and radii to create basis functions of shells
% 1 pixel thick, which are then used in the fitting.


ds=4;           % downsampling relative to raw image
cpe=16;
forceNewFile=1;
symmetrize=0;   % force a symmetrical profile
writeMiFile=1;
useAllGoodVesicles=1;  % otherwise use the vesicles selected in SimpleRSPicker
% deselectSomeVesicles=1;
includeShifts=1;
nzeros=1;
minOk=2;  % minimum value for ok field to be counted in fit.
% Modifications to the ctfs
% bScale=1;      % scaling of B-factors
% bInv=0;         % inverse B-factor
% alpha0=0.2;    %
% writeMiFile=1;

% Fitting
modelHalfwidth=32;   % angstroms
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

modelSpectrum=CCDModelSpectrum2D(1);
vm.vesicleModel=0;

tic
for fileIndex=1:numel(fname)
    disp(['Reading ' fname{fileIndex}]);
    load([infoPath fname{fileIndex}]);
    % mi.doses=mi.doses/100
    mi.basePath=rootPath;  % modify it if necessary
    if ~isfield(mi,'weights')
        mi.weights=ones(1,numel(mi.doses));
    end;
    %     iname=[mi.procPath mi.baseFilename 'm.mrc'];
    %     disp(['Reading image: ' iname]);
    m0=meReadMergedImage(mi);
    ds0=mi.imageSize(1)/size(m0,1);
    n=mi.imageSize/ds;
    m=Downsample(m0,n);
    subplot(2,2,1);
    imacs(GaussFilt(m,256/n(1)));
    drawnow;
    
    % Get the raw images
    nim=numel(mi.imageFilenames);
    fnames=cell(0);
    for i=1:nim
        fnames{i}=[mi.imagePath mi.imageFilenames{i}];
    end;
    m=meReadImages(fnames,cpe,0,mi.pixA);
    m=mePreWhiten(m,modelSpectrum)/mi.doses(1);
    [nx ny nim]=size(m);
    n0=[nx ny];
    
    %%
    if useAllGoodVesicles
        vesIndices=find(mi.vesicle.ok>=minOk);
%     if useAllVesicles
%         mi.vesicleModel=Crop(ones(27,1),33);
%         mi.vesicle.ok= mi.vesicle.r<150 & ~isnan(mi.vesicle.s);
%         vesIndices=find(mi.vesicle.ok);
%         numOk=sum(mi.vesicle.ok)
%     elseif deselectSomeVesicles
%         indices=find(mi.particle.picks(:,3)==3);  % rejected vesicles
%         disp([num2str(numel(indices)) ' vesicles not fitted']);
%         ok=ones(1,numel(mi.vesicle.x));
%         ok(mi.particle.picks(indices,4))=0;
%         vesIndices=find(ok);  % indices of vesicles to restore
    else
        disp('using only the selected ones');
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
    mi1.vesicle.s=ones(size(mi1.vesicle.s))*median(mi.vesicle.s(mi.vesicle.ok>=minOk));
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
    
    %%
    mi2=mi1;  % copy for making basis functions
    includeShifts=includeShifts && isfield(mi.vesicle,'shiftX')...
        && numel(mi.vesicle.shiftX) == numel(mi.vesicle.x);
    %     Remove oulier shifts
    if includeShifts
        nshifts=nim;
        mi2.vesicle.shiftX(abs(mi1.vesicle.shiftX)>10)=0;
        mi2.vesicle.shiftY(abs(mi1.vesicle.shiftY)>10)=0;
    else
        nshifts=1;
    end;
    R0=zeros([n nshifts nt]);   % basis without CTF
    disp(['Computing ' num2str(nt) ' x ' num2str(nshifts) ' terms']);
    
    for k=1:nshifts
        if includeShifts
            mi2.vesicle.x=mi.vesicle.x(:)+mi2.vesicle.shiftX(:,k);
            mi2.vesicle.y=mi.vesicle.y(:)+mi2.vesicle.shiftY(:,k);
        end;
        for i=1:nt
            j=(i-1-nvs)*ds+nxs+1; % pick the shell radius
            vm1=zeros(nx,1);
            vm1(j)=1;  % scale it according to thickness
            if symmetrize
                vm1(nx-j+1)=1;
            end;
            mi2.vesicleModel=vm1*ds;  % removed the mi.pixA factor
            %         mi.vesicleModel=vm1*mi.pixA;
            %         r0=meMakeModelVesicles(mi,n*ovs,vesIndices,0,0); % no filterin
            r0=meMakeModelVesicles(mi2,n,vesIndices,0,0); % no filtering
            R0(:,:,k,i)=r0;
            
            q=Crop(r0,256);
            subplot(2,2,4);
            imacs(q);
            title([i k]);
            subplot(2,2,3);
            plot(sect(q));
            drawnow;
        end;
    end;
    disp('done');
    if ~includeShifts
        R0=repmat(R0,[1 1 nim 1]);
    end;
    %%
    %     alpha0=.05;
    %     disp('filtering');
    %     for i=1:numel(mi1.ctf)
    %         mi1.ctf(i).alpha=alpha0;
    %         mi1.ctf(i).B=mi.ctf(i).B*bScale-bInv;
    %         mi.ctf(i).alpha=alpha0;  % change the original too.
    %     end;
    %
    disp('Transforming images');
    [effctf,mc,mts0,coeffz,dctfs]=meCombineImages(m,mi,ds);  % downsample
    cctf=CCDEffCTF(n,ds);
    
    %     %     Get the images ready to be summed
    %     disp('Filtering images');
    %     mtf=mts;
    %     for i=1:nim
    %         mtf(:,:,i)=real(ifftn(fftn(mts(:,:,i)).*ifftshift(coeffz(:,:,i))));
    %     end;
    %     mtf=-mtf*4;  %%%% strange.
    %     mMerge=sum(mtf,3);  % merged image for checking
    %     disp('Making model vesicles');
    %     v=meMakeModelVesicles(mi,n);
    %
    %     msub=m-v;  % basic subtracted image
    %     drawnow;
    %%
    %     figure(2);
    %     SetGrayscale;
    freqs=RadiusNorm(n)/(ds*mi.pixA);
    R=zeros(prod(n),nim,nt);  % this will be basis functions
    
    niters=50;
    niters=1;
    % Do the actual least-squares here
    vAccum=zeros(n);
    aAccum=zeros(nt,1);
    
    P=   [1  1.04  1.3    0     0    .012   .012];
    %     a1  a2    a3    al1   al2   al3     alc
    mask=[0   1      1     0     0      1      1];
    %     alc is the value inserted whenever al_n is zero.
    mi2=mi;
    
    fitWeights=[1 .5 .5];
    
    Y=zeros(prod(n),nim);
    hPlot=zeros(n(1)/2,nim);
    ctfTotal=zeros([n nim]);
    
    [coeffs mergedCTF dctfs]=meComputeMergeCoeffs( freqs, mi2.ctf, mi2.doses, nzeros, mi2.weights);
    mts=mts0;
    for i=1:nim
        q=real(ifftn(fftn(mts0(:,:,i)).*ifftshift(dctfs(:,:,i)~=0)));  % zero the data beyond zeros
        Y(:,i)=fitWeights(i)*q(:);
        mts(:,:,i)=q;
    end;
    Y2=repmat(Y,1,2);
    
    ampFactors=zeros([n nim]);
    for iter=1:niters
        if iter==niters
            P=Simplex('centroid');
        end;
        for j=1:nim  % modify ctf parameters
            %             Modify alphas
            if P(nim+j)>0
                mi2.ctf(j).alpha=P(nim+j);
            else
                mi2.ctf(j).alpha=P(2*nim+1);
            end;
            
            mi2.ctf(j).ampFactor=P(j);
            mi2.ctf(j).pixA=[];
        end;
        %         Get the overall ctfs
        [coeffs mergedCTF dctfs]=meComputeMergeCoeffs( freqs, mi2.ctf, mi2.doses, nzeros, mi2.weights);
        for j=1:nim
            %             ampFactors(:,:,j)=P(j)*(P(nim+1)+mi.ctf(j).defocus^.5/4*exp(-freqs.^2/(2*P(nim+2)^2)));
            ampFactors(:,:,j)=cctf*P(j)*fitWeights(j);
            ctfTotal(:,:,j)=ampFactors(:,:,j).*dctfs(:,:,j);
            hPlot(:,j)=sectr(ctfTotal(:,:,j))/fitWeights(j);
        end;
        subplot(4,2,6);
        plot(sectr(freqs),abs(hPlot));
        xlabel('Spatial frequency');
        ylabel('Final CTF')
        title(num2str(P));
        
        %         Filter the basis functions
        for i=1:nt
            for j=1:nim
                r=real(ifftn(fftn(R0(:,:,j,i)).*ifftshift(ctfTotal(:,:,j))));  % Filter with the ctf
                R(:,j,i)=r(:);
            end;
        end;
        R1=reshape(R,prod(n)*nim,nt);  % concatenate columns x nim.
        R1(:,nt+1)=1+0*R1(:,1);  % constant term
        a=LinLeastSquares(R1,Y);
        
        fit=R1*a;
        %         Show the result
        mFit=reshape(fit,[n nim]);
        
        d1=[0; a(1:end-1); 0];  % include zero points at the end
        d1x=meDownsampleVesicleModel(d1,1/ds);
        %     d1x=circshift(d1x,0);
        
        ndx=numel(d1x);
        if mean(d1x)<0
            d1x=-d1x;
        end;
        d1x=d1x*2/max(abs(d1x));  % force maximum to be 2.
        mi3=mi2;
        mi3.vesicleModel=d1x;
        R2=zeros([n nim 2]);
        % Now try to fit amplitudes
        for k=1:nshifts
            if includeShifts
                mi3.vesicle.x=mi.vesicle.x(:)+mi3.vesicle.shiftX(:,k);
                mi3.vesicle.y=mi.vesicle.y(:)+mi3.vesicle.shiftY(:,k);
            end;
            r0=meMakeModelVesicles(mi3,n,vesIndices,0,0); % no filtering
            for j=1:2
                mi2.ctf(k).alpha=(j-1)*pi/2;
                %         Get the overall ctfs
                [coeffs mergedCTF dctfs]=meComputeMergeCoeffs( freqs, mi2.ctf, mi2.doses, nzeros, mi2.weights);
                ctfTotal=cctf.*dctfs(:,:,k);
                
                R2(:,:,k,j)=real(ifftn(fftn(r0).*ifftshift(ctfTotal)));  % Filter with the ctf
            end;
        end;
        R3=reshape(R2,prod(n),nim*2);
        a2=LinLeastSquares(R3,Y2);
        
        fit=R3*a2;
        a2
        
        %         for i=1:nim
        %             subplot(2,2,i);
        %             imacs(GaussFilt(Crop(mts(:,:,i)-mFit(:,:,i)/fitWeights(i),1024),.1));
        %         end;
        %         subplot(4,2,8);
        %         plot([-a(1:end-1) 0*a(1:end-1)]);
        %         subplot(2,2,1);
        %         title(['Iteration ' num2str(iter)]);
        %         subplot(2,2,3);
        %         title(err);
        %         drawnow;
        %         %         pause;
        
    end;
    %%
    %     d1=a(border:end-border);  % include the zero points at the end
    d1=[0; a(1:end-1); 0];  % include zero points at the end
    if symmetrize
        d1=[d1; flipud(d1(1:nt-1))];
    end;
    d1x=meDownsampleVesicleModel(d1,1/ds);
    %     d1x=circshift(d1x,0);
    
    ndx=numel(d1x);
    if mean(d1x)<0
        d1x=-d1x;
    end;
    d1x=d1x*2/max(abs(d1x));  % force maximum to be 2.
    %     d1x=circshift(d1x,1);
    subplot(2,2,1);
    plot(1:ndx,d1x,'.-',(ndx+1)/2,d1x((ndx+1)/2),'r+','markersize',10);  % final composite profile
    xlabel([num2str(mi.pixA) ' A per pixel']);
    drawnow;
    
    %     outName=[mi.tempPath mi.baseFilename 'ds' num2str(ds) 'sym' num2str(symmetrize)...
    %         'all' num2str(useAllVesicles)...
    %         'alph' num2str(alpha0)];
    %     print('-djpeg','-r300',[outName '.jpg']);  % save the plots as a jpeg.
    %     disp('Figure image:');
    %     disp(outName);
    %
    %     vm.vesicleModel=d1x;
    %     vm.pixA=mi.pixA;
    %     save([outName '.mat'],'vm');
    if writeMiFile
        mi.vesicleModel=d1x;
        mi.ctf=mi2.ctf;
        
        disp(['Writing ' infoPath fname{fileIndex}]);
        save([infoPath fname{fileIndex}],'mi');
    end;
end;
toc
