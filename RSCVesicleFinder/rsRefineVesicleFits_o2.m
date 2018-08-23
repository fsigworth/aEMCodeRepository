% rsRefineVesicleFits2
% Revised version with acceleration options.
% Given an info structure mi, find and subtract vesicles
% and return the updated mi containing the vesicle coordinates.
% A vesicle model file is stored in Temp/ as <basename>v.mrc.
%
fittingMode=3;
% modes are:  1 - fit overall amplitude
%             2 - fit amplitude and shift only, don't change radius
%             3 - fit amplitude, shift and radius
forceNewModel=0;   % Always ask the user to select a new refined model
% on the first micrograph (can be from the same micrograph)
sRange=2;          % Outlier amplitudes must lie within this factor of the
                   % median of vesicle.s, otherwise set equal to median.
                   % Set to 1 to force median for all vesicles.
useOkField=1;      % refine every vesicle for ok(:,1) is true.
doDownsampling=1;  % Downsample for speed
maxVesiclesToFit=inf;
disA=800;          % size of displayed/fitted window in angstroms
writeMiFile=1;     % Save the updated mi file
writeSubFile=0;    % Write models and subtracted images into Temp
setBasePath=1;     % Replace the mi.basePath with the path from the file selector.
sAcceptSd=0.6;       % number of SDs to allow in ok vesicles
maxRadius=400;     % maximum radius in Å for ok to be true.
minRadius=100;
displayPeriod=5;
vm=struct;
vmGood=0;

% if nargin<1  % put up a file selector
[fname pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
if isnumeric(fname)  % Cancel
    return
end;

[rootPath infoPath]=ParsePath(pa);

if ~iscell(fname)
    fname={fname};
end;

cd(rootPath);

for fileIndex=1:numel(fname)
    %%
    disp(['Reading ' infoPath fname{fileIndex}]);
    load([infoPath fname{fileIndex}]);
    
    %     check if there is something to do.
    if numel(mi.vesicle.x)>0
        if ~isfield(mi,'log')
            mi.log=cell(0);
        end;
        
        %     Check to see if we have a generic vesicle model (central 5 points are
        %     essentially equal)
        nvm=numel(mi.vesicleModel);
        vmCtr=ceil((nvm+1)/2);
        if nvm < 3 || std(mi.vesicleModel(vmCtr-1:vmCtr+1))<1e-6...
                || forceNewModel
            if fileIndex==1  % Check the first one.
                [mname pa]=uigetfile('*mi.mat','Select an mi file for membrane model');
                if ischar(mname) % user selected something
                    disp(['Loading the membrane model from ' pa mname]);
                    vm=load([pa mname]);
                    vm=vm.mi;
                    vmGood=1;
                end;
            end;
            if vmGood  % replace the model
                disp(['Using the membrane model from ' mname]);
                if vm.pixA == mi.pixA  % Check if the same pixel size
                    mi.vesicleModel=vm.vesicleModel;  % copy the model
                else                   % resample the model
                    disp('Resampling the model');
                    mi.vesicleModel=meDownsampleVesicleModel(...
                        vm.vesicleModel,mi.pixA/vm.pixA);
                end;
            else
                disp('Using the existing generic membrane model');
            end;
        else  % the previous model was a good one.
            disp('Using the existing refined membrane model');
            vmGood=1;
            vm=mi;
        end;
               
        % Read the image and normalize to fractional contrast
        [m mergePath]=meReadMergedImage(mi);
        % Estimate the high-frequency variance
        sp1=mean(RadialPowerSpectrum(m),2);
        nsp=numel(sp1);
        hfVar0=mean(sp1(round(0.6*nsp):round(0.95*nsp)));
        
        %     Check that we have a temp directory
        if ~isfield(mi,'tempPath');
            mi.tempPath='Temp/';
        end;
        if ~exist(mi.tempPath,'dir')
            mkdir('Temp');
        end;
        if setBasePath
            mi.basePath=rootPath;
            disp(['Resetting mi.basePath to: ' mi.basePath]);
        end;

        mi1=mi;  % save the original in mi, and load the new model into mi1
        
        % Get image and pixel sizes
        n=size(m,1);
        ds0=mi.imageSize(1)/n;  % downsampling factor of m
        pixA0=mi.pixA*ds0;    % pixel size of m
        if doDownsampling
            % downsample the merged image to about 10A per pixel, yielding the image ms
            targetPixA=12;  % maximum pixel size
            ns=NextNiceNumber(n*pixA0/targetPixA,5,4);  % multiple of 4, max factor 5.
            if ns<n
                disp(['Downsampling to ' num2str(ns) ' pixels.']);
                ms=Downsample(m,ns);
            else
                ns=n;
                ms=m;
            end;
            ds=ds0*n/ns;  % downsampling factor of ms relative to original images.
            pixA=ds*mi.pixA;  % pixA in the image ms.
        else  % use the original merged image scale
            ds=ds0;
            pixA=pixA0;
            ns=n;
            ms=m;
        end;
        hfVar=hfVar0*(ds0/ds)^2;  % hf spectral density after downsampling
        ndis=NextNiceNumber(disA/pixA);  % size of display/fitting image
        
        %%  Get the original subtraction, and modify the amplitudes if necessary.\
        % subtract every vesicle that was originally fit
        mi1.vesicle.s(isnan(mi.vesicle.s))=0;
        sMedian=median(mi1.vesicle.s);
        %   Force ridiculous amplitude values to the median.
        q=mi1.vesicle.s<sMedian/sRange | mi1.vesicle.s > sMedian*sRange;
        mi1.vesicle.s(q)=sMedian;
        %   Compute the basic (old) subtraction using mi1 with new model
        vs=meMakeModelVesicles(mi1,ns,find(mi.vesicle.ok(:,1)));
        %   Use linear least-squares to fine-tune the scaling
        msAmpScale=(vs(:)'*ms(:))/(vs(:)'*vs(:))
        if msAmpScale>1e-3 % don't allow ridiculous values
            vs=vs*msAmpScale;
            mi1.vesicle.s=mi1.vesicle.s*msAmpScale;
        end;
        msub=ms-vs;  % This is the subtracted image we'll use

        %         Whiten and mask the subtracted image
        mmask=meGetMask(mi1,ns);
        pwH=meGetNoiseWhiteningFilter(mi1,ns);
        msubf=mmask.*real(ifftn(fftn(msub).*ifftshift(pwH)));
        figure(1);
        SetGrayscale;
        subplot(2,3,1);
        imacs(msubf);
        title(['Old subtraction, scaling ' num2str(msAmpScale)]);
        
        %%
        switch fittingMode
            case 1  % normalize the amplitude of the whole image only
%                 vs1=meMakeModelVesicles(mi1,ns).*mmask;  % new subtraction
%                 ms1AmpScale=(vs1(:)'*ms(:))/(vs1(:)'*vs1(:))
%                 if ms1AmpScale>1e-3 % don't allow ridiculous values
%                     vs1=vs1*ms1AmpScale;
%                 end;
%                 msub1=ms-vs1;
                subplot(2,3,2);
                imacs(msub);
                title('New subtraction');
%                 mi1.vesicle.s=mi1.vesicle.s*ms1AmpScale;
                q=isnan(mi.vesicle.s);
                mi1.vesicle.s(q)=0;
                mi1.vesicle.ok(:,3)=mi1.vesicle.ok(:,1) & ~q;

            case {2 3}  % tune up the x, y and s for each vesicle.
                %
                %                 if ~isfield(mi.vesicle,'ok')
                %                     mi1.vesicle.ok=ones(size(mi1.vesicle.x,3));
                %                 end;
                %%  Actual fitting is done here
                vfit=zeros(ns,ns);
                figure(1)
                disp('Fine fitting translation and amplitude.');
                drawnow;
                nVesicles=sum(mi1.vesicle.ok(:,1)>0);
                nVesicles=min(nVesicles,maxVesiclesToFit);
                nVesicles
                effCTF=ifftshift(meGetEffectiveCTF(mi1,ndis,ds));
                pwFilter=ifftshift(meGetNoiseWhiteningFilter(mi1,ndis,ds));
                figure(1);
                mi1.vesicle.ok(:,3)=true;  % we'll mark unfitted vesicles here.
                for ind=1:nVesicles
                    ok=mi1.vesicle.ok(ind,1)>0;  % The vesicle exists
                    if ~useOkField || ok
                        doDisplay=mod(ind,displayPeriod)==0;
                        [mi1 diffIm vesFit]=rsQuickFitVesicle(msubf,mmask,mi1,mi1,...
                            ind,effCTF.*pwFilter,hfVar,fittingMode,doDisplay);
                        
                        %                         subplot(2,2,2);
                        %                         title(num2str([ind mi1.vesicle.s(ind) mi.vesicle.ok(ind,2)]));
                        %                         pause;
                        %                         title(fname{fileIndex},'interpreter','none');
                    else
                        mi1.vesicle.s(ind)=NaN;
                    end;
                end;
                maxR=maxRadius/mi1.pixA;
                minR=minRadius/mi1.pixA;
                q=isnan(mi1.vesicle.s);  % blank the unfittable vesicles.
                mi1.vesicle.s(q)=0;
                mi1.vesicle.ok(:,3)=~q;  % unrefinable vesicles are marked 0
                %                 medianS=median(mi1.vesicle.s);
                %                 good=(mi1.vesicle.s>medianS/(1+sAcceptSd)) & (mi1.vesicle.s<(1+sAcceptSd)*medianS)...
                %                     & (mi1.vesicle.r < maxR) & (mi1.vesicle.r > minR);
                %              oldOk=mi.vesicle.ok;
                %                 mi1.vesicle.ok=min(mi1.vesicle.ok,1+good);
                numberRefined=sum(mi1.vesicle.ok(:,3))
                numberGood=sum(all(mi1.vesicle.ok(:,1:3),2))  %    exists, in range, refined
                mi1.vesicle.refined=1;
                
        end; % switch
        
        %%  Outputting
        mi=mi1;
        if writeMiFile
            %%
            outName=[infoPath mi.baseFilename 'mi.mat'];  % replace original
            if ~isfield(mi,'log')
                mi.log=cell(0);
            end;
            mi.log{numel(mi.log)+1}=['meRefineVesicleFits ' TimeStamp];
            save(outName,'mi');
            disp([outName ' saved']);
            
        end;
        %%             Compute and store model vesicles
        figure(2); clf; SetGrayscale;
        if writeSubFile
            
            imacs(GaussFilt(m,.1));
            title('Original image');
            drawnow;
            
            disp('Making the final vesicle models');
            vs1=meMakeModelVesicles(mi,n,find(mi.vesicle.ok(:,3)));
            
            imacs(GaussFilt(m-vs1,.1));
            title('Subtracted');
            drawnow;
            outVesName=[mi.tempPath mi.baseFilename 'v'];
            WriteMRC(vs1,pixA0,[outVesName '.mrc']);
            WriteJpeg(vs1,outVesName);
            imwrite(uint8(imscale(rot90(vs1),256,0)),[outVesName '.jpg']);
            disp([outVesName ' saved']);
            %               Store the subtracted micrograph
            outSubName=[mi.tempPath mi.baseFilename 'mv'];
            WriteMRC(m-vs1,pixA0,[outSubName '.mrc']);
            WriteJpeg(m-vs1,outSubName);
            %             imwrite(uint8(imscale(rot90(m-vs1),256,1e-3)),[outSubName '.jpg']);
            disp([outSubName ' saved']);
        else
            imacs(ms);
            title('Original');
            drawnow;
            disp('Making the final vesicle models');
            vs1=meMakeModelVesicles(mi,n,find(mi.vesicle.ok(:,1)));
            imacs(ms-vs);
            title('Subtracted');
            drawnow;
        end;
        
    else  % No vesicles have been found to refine
        if numel(mi.vesicle.x)<1
            disp('  ...no vesicles found');
        end;
        if isfield(mi.vesicle,'refined') && mi.vesicle.refined>0
            disp('  ...already refined.');
        end;
    end;
    disp(' ');
end; % for fileIndex
disp('Done.');
disp(' ');
