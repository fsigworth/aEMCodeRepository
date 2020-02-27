% rsRefineMembraneModel2
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
useAllGoodVesicles=1;  % otherwise only use the vesicles selected in SimpleRSPicker
includeShifts=1;
nzeros=1;
% Modifications to the ctfs
% bScale=1;      % scaling of B-factors
% bInv=0;         % inverse B-factor
% alpha0=0.2;    %
% writeMiFile=1;

% Fitting
modelHalfwidth=32;   % angstroms
modelHalfwidth=36;   % angstroms
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

tic
for fileIndex=1:nfiles
    disp(['Reading ' fname{fileIndex}]);
    load([infoPath fname{fileIndex}]);
    mi.basePath=rootPath;  % modify it if necessary
    if ~isfield(mi,'weights')
        mi.weights=ones(1,numel(mi.doses));
    end;
    
    % Get the data images from the original micrograph files
    mts0=meMakeMergeImageSet(mi,cpe,ds);
    [nx,ny,nim]=size(mts0);
    n=[nx ny];

%     Display the merged image
    subplot(2,2,1);
    imacs(GaussFilt(sum(mts0,3),256/nx));
    drawnow;
    
    %%
    if useAllGoodVesicles
        vesIndices=find(all(mi.vesicle.ok,2));
    else
        disp('using only the selected ones');
        vesPicks=find(mi.particle.picks(:,3)==2);  % type=2 means a picked vesicle.
        vesIndices=mi.particle.picks(vesPicks,4);
        nv=numel(vesIndices);
        disp([num2str(nv) ' vesicles to be fitted']);
    end;
    
    % Create an mi structure for modification
    mi1=mi;
    %         Force all the vesicles to the same amplitude
    mi1.vesicle.s=ones(size(mi1.vesicle.s))*median(mi.vesicle.s(vesIndices));

    % Decide how many terms to compute
    nvs=round(modelHalfwidth/(ds*mi.pixA));
    nxs=nvs*ds;  % halfwidth of the oversampled model.
    nm=2*nvs+1;  % number of points in the model
    nx=2*nxs+1;
    if symmetrize
        nt=nvs+1;  % number of terms in the model fit
    else
        nt=nm;
    end;
    
    %%  Compute the nt x nim images of shells
    %    R0(x,y,nim,nt)
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
        vm1=zeros(nx,nt);
        for i=1:nt
            j=(i-1-nvs)*ds+nxs+1; % pick the shell radius
            vm1=zeros(nx,1);
            vm1(j,i)=1;
            mi2.vesicleModel=vm1(:,i)*ds;  % removed the mi.pixA factor
            r0=meMakeModelVesicles(mi2,n,vesIndices,0,0); % no filtering
            R0(:,:,k,i)=r0;
            
%             Show an enlarged view of the center of the shell models
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

    %%  Do the big least-squares fit

    freqs=RadiusNorm(n)/(ds*mi.pixA);
    R=zeros(prod(n),nim,nt);  % this will be basis functions for ls fit
    ccdTransferFunction=CCDEffCTF(n,ds);
    
    niters=50;

%     vAccum=zeros(n);
%     aAccum=zeros(nt,1);
    
    P=   [1  1.04  1.3    0     0    .012   .012  100  100];
    %     a1  a2    a3    al1   al2   al3     alc   B0  Bd
    mask=[0   1      1    0     0      1      1    0    0 ];
    modifyBs=0;
    %     alc is the value inserted whenever al_n is zero.
    mi2=mi;
    for i=1:nim
        mi2.ctf(i).pixA=[];  % get rid of this field
    end;
    fitWeights=[1 .15 .15];  % relative weight to give each image in ls fit
    
    Y=zeros(prod(n),nim);
    hPlot=zeros(n(1)/2,nim);
    ctfTotal=zeros([n nim]);
    bs=zeros(2,nim);
    const=ones(prod(n),1);
    
%   Get the radiation-decay weighted dctfs, no ampFactor
    [coeffs, mergedCTF, dctfs]=meComputeMergeCoeffs( freqs, mi2, nzeros);
    mts=mts0;
    for i=1:nim % make Y in the ls problem Ra=Y.  It is the concatenation
%                 of the nim images, after zeroing them 
%                 beyond the nth zero, but don't phase-flip.
%         q=real(ifftn(fftn(mts0(:,:,i)).*ifftshift(coeffs(:,:,i)~=0)));  % zero the data beyond zeros
        q=mts0(:,:,i);
        Y(:,i)=fitWeights(i)*q(:);
        mts(:,:,i)=q;
    end;
    
    P=Simplex('init',P,0.2,mask);
%
    for iter=1:niters
        if iter==niters
            P=Simplex('centroid');
        end;
        for j=1:nim  % modify ctf parameters
            if modifyBs
                mi2.ctf(j).B=P(2*nim+2)+mi.ctf(j).defocus*P(2*nim+3);
            end;
            %             Modify alphas
            if P(nim+j)>0
                mi2.ctf(j).alpha=P(nim+j);    % Explicit value given
            else
                mi2.ctf(j).alpha=P(2*nim+1);  % default for all others
            end;
            
            mi2.ctf(j).ampFactor=P(j);  % set the amplitude of this exposure
        end;
        %         Get the ctf for each image (theoretical ctf * rad damage)
        [coeffs, mergedCTF, dctfs]=meComputeMergeCoeffs( freqs, mi2.ctf, mi2.doses, nzeros, mi2.weights);
        for j=1:nim
            %  Get the overall ctf for each image, including ccd and ampFactors 
            ctfTotal(:,:,j)=-dctfs(:,:,j).*coeffs(:,:,j).*ccdTransferFunction*mi2.ctf(j).ampFactor*fitWeights(j);
            hPlot(:,j)=sectr(ctfTotal(:,:,j))/fitWeights(j);
        end;
%             Plot the overall ctfs
        subplot(4,2,6);
        plot(sectr(freqs),abs(hPlot));
        xlabel('Spatial frequency');
        ylabel('Final CTF')
        title(num2str(P));
        
        %         Filter the basis functions according to the overall ctf
        for i=1:nt  % loop over terms
            for j=1:nim  % loop over images
                r=real(ifftn(fftn(R0(:,:,j,i)).*ifftshift(ctfTotal(:,:,j))));  % Filter with the ctf
                R(:,j,i)=r(:);
            end;
        end;
         R1=reshape(R,prod(n)*nim,nt);  % concatenate columns x nim.
         R1(:,nt+1)=1+0*R1(:,1);  % append a constant term
%         least-squares fit of coefficients
        a=LinLeastSquares(R1,Y);
%         Simplex optimization of the ctf parameters
        fit=R1*a;
        diff=(fit(:)-Y(:));
        err=diff'*diff;
        err=err+1e6*sum(P<0);  % don't allow negative values.
        if mod(iter,10)==0
            disp([P err/1e5]);
        end;
        P=Simplex(err);
        
        %         Show the results
        mFit=reshape(fit,[n nim]);
        for i=1:nim
            subplot(2,2,i);
            imacs(GaussFilt(mts(:,:,i)-mFit(:,:,i)/fitWeights(i),.1));
        end;
    
        subplot(4,2,8);
        plot([-a(1:end-1) 0*a(1:end-1)]);
        subplot(2,2,1);
        title(['Iteration ' num2str(iter)]);
        subplot(2,2,3);
        title(err);
        drawnow;
        
    end;
    %%  Put the membrane profile into standard form
    %     d1=a(border:end-border);  % include the zero points at the end
%     d1=[0; a(1:end-1); 0];  % include zero points at the end
    d1=[0; a; 0];  % include zero points at the end; not using const.
    if symmetrize
        d1=[d1; flipud(d1(1:nt-1))];
    end;
    d1x=meDownsampleVesicleModel(d1,1/ds);
    
    ndx=numel(d1x);
    if mean(d1x)<0
        d1x=-d1x;
    end;
    d1x=d1x*2/max(abs(d1x));  % force maximum to be 2.

    subplot(2,2,1);
    plot(1:ndx,d1x,'.-',(ndx+1)/2,d1x((ndx+1)/2),'r+','markersize',10);  % final composite profile
    xlabel([num2str(mi.pixA) ' A per pixel']);
    drawnow;
    
    if writeMiFile
        mi.vesicleModel=d1x;
        mi.ctf=mi2.ctf;
        
        disp(['Writing ' infoPath fname{fileIndex}]);
        save([infoPath fname{fileIndex}],'mi');
    end;
end;
toc
