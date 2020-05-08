% rsRefineMembraneModel4
% Fit the differently-defocused images separately.
% Fit multiple micrographs simultaneously.
%
% Use linear least-squares fitting to deduce a membrane profile from the
% vesicles in a group of images.  Given an info structure mi and a merged image, use
% the stored vesicle centers and radii to create basis functions of shells
% 1 pixel thick, which are then used in the fitting.
% Write out the identical membrane model in each mi file of the group.
% Repeat this for multiple groups.

groupSize=28;    % make this large to average all the fits
% ds=2;           % downsampling relative to raw image
ds=4;           % downsampling relative to raw image
% ds=6;           % downsampling relative to raw image
nIters=20;
fitWeights=[1 .5 .4];  % relative weight to give each image in ls fit
fHP=.003;   % Highpass filter freq. in A^-1

vesMinR=150;  % in Å
vesMaxR=300;
vesAmpFactor=1.5;  % ratio of s(:,1) relative to median
vesMaxEllip=3;  % ellipsitiy in Å
vesMaxShadow=.2; % rel density

% for 2 exposures the middle element is deleted.

% *** Special values to fit mainly the 1st exposure ***
% nIters=20;
% fitWeights=[1 .3 .1];  % %%%% relative weight to give each image in ls fit
addP4=0;  % special insertion of 'p4' into imageFilenames
forceNewFile=0; % always put up a file selector, don't repeat the previous files
forceMedianAmplitude=0;  % force vesicle.s to be set to global median
symmetrize=0;   % force a symmetrical profile
writeMiFile=1;  % Write out the results
defaultCpe=16;
maxShift=10;
valueB=20;


useAllGoodVesicles=1;  % otherwise only use the vesicles selected in SimpleRSPicker
includeShifts=0; % Use data from meTrackVesicleShifts
nzeros=1;       % Use high-defocus data up to this number of zeros

% Fitting
% modelHalfwidth=36;   % angstroms
modelHalfwidth=30;   % angstroms
%  modelHalfwidth=42;
fitDC=1;        % also fit a constant offset to the image.

figure(1); clf; SetGrayscale;

% If we've already run this, we don't ask for a filename
if ~(exist('allNames','var') && exist('rootPath','var')) || forceNewFile  % not already defined
    [allNames, pa]=uigetfile('*mi.*','Select mi files','multiselect','on');
    if isnumeric(allNames)
        return
    end;
    % pa
    [rootPath, infoPath]=ParsePath(pa);
    if ~iscell(allNames)
        allNames={allNames};
    end;
end;
cd(rootPath);
totalFiles=numel(allNames);

% We will fit the files in groups
ngroups=max(1,floor(totalFiles/groupSize))

tic

%%  Huge loop over groups

for groupIndex=1:ngroups
    mis=cell(0);  % copy of all the mi structures
    mis1=cell(0); % modified mi structures for fits
    startFileIndex=1+(groupIndex-1)*groupSize;
    endFileIndex=min(groupIndex*groupSize,totalFiles);
    if groupIndex==ngroups
        endFileIndex=totalFiles;
    end;
    disp(' ');
    disp(['Group ' num2str(groupIndex) '  Files: ' num2str([startFileIndex endFileIndex])]);
    
    fname = allNames(startFileIndex:endFileIndex);
    nfiles=endFileIndex-startFileIndex+1;
    
    for h=1:nfiles
        disp(['Reading ' fname{h}]);
        mi=ReadMiFile([infoPath fname{h}],1);
        mi.basePath=rootPath;  % modify it if necessary
        if ~isfield(mi,'weights')
            mi.weights=ones(1,numel(mi.doses));
        end;
        if ~isfield(mi,'cpe')
            mi.cpe=defaultCpe;
            disp(['mi.cpe set to ' num2str(mi.cpe)]);
        end;
        if addP4
            for i=1:numel(mi.imageFilenames)
                str=mi.imageFilenames{i};
                if ~exist([mi.imagePath str],'file')
                    q=strfind(str,'sum');
                    if numel(q)>0
                        q=q(1);
                        mi.imageFilenames{i}=[str(1:q-1) 'p4' str(q:end)];
                    end;
                end;
            end;
        end;
        mis{h}=mi;
        % Get the data images from the original micrograph files
        mts=meMakeMergeImageSet(mi,mi.cpe,ds);
        if h==1  % initialize
            [nx,ny,nim]=size(mts);
            n=[nx ny];
            % Decide how many terms to compute
            pixA=mi.pixA*ds;  % working pixel size
        hpKernel=GaussHPKernel(n,fHP*pixA);
    freqs=RadiusNorm(n)/pixA;
        ccdTransferFunction=CCDEffCTF(mi.camera,n);
        
            nvs=round(modelHalfwidth/pixA);
            nxs=nvs*ds;  % halfwidth of the oversampled model.
            nm=2*nvs+1;  % number of points in the model
            nx=2*nxs+1;
            if symmetrize
                nt=nvs+1;  % number of terms in the model fit
            else
                nt=nm;
            end;
            mtsi=zeros([n nim nfiles],'single');  % array of all the images
            maski=false([n nfiles]);  % Mask for each merged image
            Ri=zeros([n nfiles nim nt],'single');  % big array for basis functions
            %         R=zeros([n nfiles nim nt],'single');  % big reference array
        end;

        
        %%         Collect all the images and masks
        mtsHP=GaussHP(mts,fHP*pixA,1);
        mtsi(:,:,:,h)=mtsHP;
        maski(:,:,h)=meGetMask(mi,n);
        %     Display the merged images
        subplot(2,2,1);
        imacs(GaussFilt(sum(mts,3),256/nx));
        title(fname{h},'interpreter','none');
        drawnow;
        
        %
        if isfield(mi.vesicle,'x') && numel(mi.vesicle.x)>0 % there are some
        if useAllGoodVesicles
            vesOk=(all(mi.vesicle.ok,2)); % marked good, and refined
%           check that the amplitude is near the median
            medS=median(mi.vesicle.s(vesOk,1));
            vesOk=vesOk & mi.vesicle.s(:,1)>medS/vesAmpFactor & mi.vesicle.s(:,1)<medS*vesAmpFactor;
%             Check the geometric criteria
            if size(mi.vesicle.r,2)<3
                mi.vesicle.r(1,3)=0;
            end;
            if size(mi.vesicle.s,2)<2
                mi.vesicle.s(1,2)=0;
            end;
            vesOk=vesOk & (mi.vesicle.r(:,1)>vesMinR/mi.pixA...
                & mi.vesicle.r(:,1)<vesMaxR/mi.pixA...
                & abs(mi.vesicle.r(:,3))<vesMaxEllip/mi.pixA...
                & abs(mi.vesicle.s(:,2))<vesMaxShadow*mi.vesicle.s(:,1));
%             disp([num2str(sum(vesOk)) ' vesicles out of ' num2str(numel(vesOk)) ' used.']);
            vesIndices=find(vesOk);
        else
            disp('Using only the selected vesicles');
            vesPicks=find(mi.particle.picks(:,3)==2);  % type=2 means a picked vesicle.
            vesIndices=mi.particle.picks(vesPicks,4);
        end;
        else
            vesIndices=[];
        end;
            nv=numel(vesIndices);
            disp([num2str(nv) ' vesicles to be fitted out of ' num2str(numel(mi.vesicle.x))]);

        % Create an mi structure for modification
        mi1=mi;
        mi1.vesicle.extraPeaks=[];
        %         Force all the vesicles to the same amplitude
        if forceMedianAmplitude
            mi1.vesicle.s=ones(size(mi1.vesicle.s))*median(mi.vesicle.s(vesIndices));
        end;
        
        %%  Compute the nim x nfiles x nt images of shells
        %    R0(x,y,nim,nfiles nt)
        for i=1:numel(mi1.ctf)
            mi1.ctf(i).ampFactor=1;
        end;
        mi2=mi1;  % copy for making basis functions
        includeShifts=includeShifts && isfield(mi.vesicle,'shiftX')...
            && numel(mi.vesicle.shiftX) == numel(mi.vesicle.x);
        %     Remove oulier shifts
        if includeShifts
            nshifts=nim;
            mi2.vesicle.shiftX(abs(mi1.vesicle.shiftX)>maxShift)=0;
            mi2.vesicle.shiftY(abs(mi1.vesicle.shiftY)>maxShift)=0;
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
            for i=1:nt  % loop over shells
                j=(i-1-nvs)*ds+nxs+1; % pick the shell radius
                vm1=zeros(nx,1);
                vm1(j)=1;
                mi2.vesicleModel=vm1*ds;  % removed the mi.pixA factor
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
        if ~includeShifts
            R0=repmat(R0,[1 1 nim 1]);
        end;
        
        Ri(:,:,h,:,:)=R0;
    end; % loop over files indexed by h.
    %%  Do the big least-squares fit
    
    R=zeros(prod(n),nfiles,nim,nt);  % this will be basis functions for ls fit
    
    switch nim
        case 3
            P=   [1  1.04  1.3    0     0      0   .02  50  50];
            %     a1  a2    a3    al1   al2   al3   alc   B0  Bd
            mask=[0   0      0     0     0     0     1    0    0 ];
        case 2
            %             P=   [1  1.3     0     0    .01  100  100];
%             P=   [1  1.3     0     0    .01  50  50];
            P=   [1  1.3     0     0    .01  valueB valueB];
            %     a1  a2    al1   al2   alc   B0  Bd
            mask=[0   0      0     0     1    0    0 ];
        case 1
            P=   [1  1.3     0     0    .01  valueB valueB];
            %     a1  a2    al1   al2   alc   B0  Bd
            mask=[0   0      0     0     1    0    0 ];
otherwise
            warning('Fit params inconsistent with number of exposures');
    end;
    while numel(fitWeights)>nim
        fitWeights(2)=[];
    end;
    %     alc is the value inserted whenever al_n is zero.
    modifyBs=any(P(2*nim+2:2*nim+3))~=0;
    endNullFactor=3e6;
    Y=zeros(prod(n),nfiles,nim);
    hPlot=zeros(n(1)/2,nim);
    ctfTotal=zeros([n nim]);
    bs=zeros(2,nim);
    const=ones(prod(n)*nfiles,1);
    
    % %   Get the radiation-decay weighted dctfs, no ampFactor
    %     [coeffs, mergedCTF, dctfs]=meComputeMergeCoeffs( freqs, mi2, nzeros);
    for h=1:nfiles
        mts=mtsi(:,:,:,h);
        msk=maski(:,:,h);
        for i=1:nim % make Y in the ls problem Ra=Y.  It is the concatenation
            %                 of the nim images
            q=mts(:,:,i).*msk;  % it is masked like the basis fcns
            Y(:,h,i)=fitWeights(i)*q(:);
        end;
    end;
    Ps=P;
    P=Simplex('init',P,0.7,mask);
    pOld=P;
    for iter=1:nIters
        Ps(iter,:)=P;
        
        for h=1:nfiles
            mi2=mis{h};
            for j=1:nim  % modify ctf parameters
                if modifyBs  % otherwise, leave ctf B values untouched
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
            [coeffs, mergedCTF, dctfs]=meComputeMergeCoeffs2(...
                freqs, mi2, nzeros);
            for j=1:nim
                %  Get the overall ctf for each image, including ccd and ampFactors
                ctfTotal(:,:,j)=hpKernel.*coeffs(:,:,j).*dctfs(:,:,j).*ccdTransferFunction...
                    *mi2.ctf(j).ampFactor*fitWeights(j)*mi2.weights(j);
                hPlot(:,j)=sectr(ctfTotal(:,:,j))/fitWeights(j);
            end;
            %             Plot the overall ctfs
            subplot(4,2,6);
            plot(sectr(freqs),abs(hPlot));
            xlabel('Spatial frequency');
            ylabel('Final CTF')
            title(num2str([h P]));
            drawnow;
            msk=maski(:,:,h);
            %         Filter the basis functions according to the overall ctf
            %         and mask each one.
            for i=1:nt  % loop over terms
                for j=1:nim  % loop over images
                    %                     Filter each refernce with the ctf, and mask it.
                    r=real(ifftn( fHP.*fftn(Ri(:,:,h,j,i)).*ifftshift(ctfTotal(:,:,j)) )).*msk;
                    R(:,h,j,i)=r(:);
                end;
            end;
        end; % for h
        R1=double(reshape(R,prod(n)*nfiles*nim,nt));  % concatenate columns x nim.
        %         Y is [n x nfiles x nim.  both Y and R1 are masked.
        %         Least-squares fit to get the nt membrane profiles
        a=LinLeastSquares(R1,Y(:));
        fit=single(R1*a);  % [n x nfiles x nim]
        
        %         Show the results
        mFit=reshape(fit,[n nfiles nim]);
        for i=1:nim
            subplot(2,2,i);  % display the parts of the first micrograph
            m=reshape(mtsi(:,:,i,1),n);
            imacs(GaussFilt((m-mFit(:,:,1,i)/fitWeights(i)).*maski(:,:,1),.1));
            title(num2str([h i]));
            drawnow;
        end;
        
        %   Do a separate ls fit for amplitude factor for each of the nim
        %   exposures.
        models=reshape(fit,prod(n)*nfiles,nim);
        err=0;
        Y1=reshape(Y,prod(n)*nfiles,nim);
        for i=1:nim
            F=[models(:,i) const];  % include a constant offset
            y=Y1(:,i);
            b=LinLeastSquares(double(F),double(y));
            bs(:,i)=b;
            diff=(F*b-y);
            err=err+diff'*diff;
            Y(:,i)=Y(:,i)-bs(2,i);  % correct the baseline
        end;
        P(1:nim)=P(1:nim).*bs(1,:)/bs(1,1);  % insert the amplitude parameters
        diff=(fit(:)-Y(:));
        %         err=err+endNullFactor*sum(abs([a(1:2)' a(end-1:end)']));
        err=err+1e6*sum(P<0);  % don't allow negative values.
        %     if mod(iter,10)==0
        %         disp([P err/1e5]);
        %     end;
        pOld=P;  % we'll reserve the amplitude values.
        if iter==nIters  % last iterations
            P=Simplex('centroid');
        else
            P=Simplex(err);
        end;
        P(1:nim)=pOld(1:nim);  % keep the old amplitudes
        
        subplot(4,2,8);
        %         plot([-a(1:end-1) 0*a(1:end-1)]);
        plot([a 0*a]);
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
    d1x=meDownsampleVesicleModel(d1,1/ds);
    
    ndx=numel(d1x);
    if mean(d1x)<0
        d1x=-d1x;
    end;
    d1x=d1x*2/max(abs(d1x));  % force maximum to be 2.
    
    subplot(2,2,3);
    plot(1:ndx,d1x,'.-',(ndx+1)/2,d1x((ndx+1)/2),'r+','markersize',10);  % final composite profile
    xlabel([num2str(mi.pixA) ' A per pixel']);
    title(fname{1},'interpreter','none');
    drawnow;
    %%
    %         Save the figure as a jpeg
    tempDir='Temp/';
    if ~exist(tempDir,'dir')
        mkdir(tempDir);
    end;
    set(gcf,'paperpositionmode','auto');
    jName=[tempDir 'MbnG' num2str(groupSize) mi.baseFilename '.jpg'];
    print('-djpeg','-r300',jName);
    %%
    if writeMiFile
        disp('Saving files:');
        
        for h=1:nfiles
            mi2=mis{h};
            for j=1:nim  % modify ctf parameters
                if modifyBs  % otherwise, leave ctf B values untouched
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
            mi2.vesicleModel=d1x;
            
            mi=mi2;
            mi.log{end+1,1}=['rsRefineMembraneModel4 ' TimeStamp];

            disp([infoPath fname{h}]);
            WriteMiFile(mi,[infoPath fname{h}]);
        end;
        
        
        %         disp(['Writing ' infoPath fname{h}]);
        %         save([infoPath fname{h}],'mi');
    end;
end; % for groupIndex
toc

% Show some vesicle statistics
figure(2);


