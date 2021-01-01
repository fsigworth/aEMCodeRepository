function rsRefineVesicleFits(miNames,mpars)
% function rsRefineVesicleFits(miNames,mpars)
% Fits distorted vesicles.
% Revised version with acceleration options.
% miNames is a cell array of mi file names (with path)
% mpars is a struct of parameters (see the default values set below)
% Find and subtract vesicles
% and return the updated mi containing the vesicle coordinates, radius
% expansion and amplitude angular expansion.
% If desired, a vesicle model file is stored in Temp/ as <basename>v.mrc.
% If desired, a subtracted file is stored in Merged/ as <basename>mv.mrc.

% if nargin<1
%     miNames=[];
% end;
disp('rsRefineVesicleFits startup.');


if nargin<2
    mpars=struct;
end;
% --------Parameter defaults--------
pa=AddSlash(fileparts(which('rsRefineVesicleFits'))); % Get our directory
dpars.modelMiName=[pa 'ModelMi.txt'];

dpars.overwrite=1;
dpars.writeMiFile=1;
dpars.doPreSubtraction=1; % subtract all vesicles except the one being fitted
dpars.listFits=1;  % print out each fit's parameters
dpars.scaleOriginalAmplitudes=1; % multiply amplitudes by this value

dpars.writeSubMRC=1; % Write out subtracted image <basename>mv.mrc

dpars.writeSmallMRC=0; % Write out image downsampled to M4. Otherwise, we 
% match the size of any existing small image e.g. Merged_sm/*ms.mrc when we
% write a small sub mrc.
dpars.writeSmallSubMRC=1; % Write out a downsampled subtracted image
dpars.dsSmall=4; % downsampling factor for small output images
% dpars.writeVesFiles=0;   % Write vesicle models into Temp/

dpars.resetBasePath=1;   % update the basePath field of the mi file to the current directory.
dpars.modifiedAlpha=0;  % Good alpha value for lipids.

% Number of terms (for both radius and amplitude fitting) is set thusly:
%    nTerms=find(vesicle.r(ind,1) < pars.rTerms/mi.pixA,1);
% i.e. nTerms is the index of last entry of rTerms smaller than our radius.
% dpars.rTerms=[100 150 200 250 250 300 inf];
dpars.rTerms=[100 150 200 300 inf];  % quick
% dpars.rTerms=[150 200 inf];  % for Mengqiu

% Define the fitting mode for each round
% dpars.fitModes={'RadiusOnly' 'RadiusAndLin'};
dpars.fitModes={'RadiusOnly' 'LinOnly'}; % First fit shape, then amplitudes
% dpars.fitModes={'LinOnly'};
%    dpars.fitModes={'RadiusAndLin'};
dpars.fractionStartingTerms=[.5 1]; % in radius-fitting modes, 
%    thefraction of terms to use in each round
dpars.fractionAmpTerms=[0 1]; % fraction of amp terms to use
%  To avoid crazy fits, we repeat the whole radius-only fitting with the base
%  radius perturbed by these steps (in angstroms) and pick the best.
dpars.radiusStepsA=[-100 -50 0 50]; % repeat radius-only fitting with perturbed r(1)
dpars.maxPixA=6;  % downsampled image resolution for radius fitting
dpars.maxPixA=3; %%%%%%%%%%%

dpars.disA=1200;  % size of the fit window, in angstroms
% dpars.disA=1600;  % for Mengqiu

%  We add extra peaks in the scattering profile when fitting amplitudes
dpars.peakPositionA=[-37 0 37];  % empirical default.  Works a bit better than [37 37]
dpars.xPeakSigmaA={5 5}; % width of extra Gaussian peaks, in angstrom

% -----------Merge the defaults with the given mpars-----------
pars=SetOptionValues(dpars,mpars,1);

%     The following must have at least as many elements as dpars.fitModes!
pars.xPeakPositionA=cell(1,numel(pars.fitModes));
pars.xPeakPositionA{end}=pars.peakPositionA; % add the peaks to the last round.

% -----------------------

% See if we are on the cluster.
host=getenv('HOSTNAME');
disp(['host name: ' host]);
isCluster=strncmp(host,'c',1); % if we are on a farnam node
batchMode=isCluster  && (exist('miNames','var') && numel(miNames)>0);

displayOn=~batchMode;
% displayOn=0; % 1: show results at end of each fit; 2 update while fitting.

dsModel=2;  % net downsampling of "full res" model relative to orig micrograph
dsSmall=pars.dsSmall;  % net downsampling of 'small' output file.

forceNewModel=0;   % Always ask the user to select a new refined model
  % on the first micrograph (can be from the same micrograph)
resetBasePath=pars.resetBasePath;
writeMiFile=pars.writeMiFile;     % Save the updated mi file

vm=struct;         % our local vesicle (membrane) model
vmGood=0;          % model is valid

% If none given, have the user select some mi files: boilerplate
if ~(exist('miNames','var') && numel(miNames)>0)
    [miNames, pa]=uigetfile({'*mi.txt' '*mi.mat'},'Select mi files','multiselect','on');
    if isnumeric(pa) % File selection cancelled
        return
    end;
    [rootPath, infoPath]=ParsePath(pa);
    cd(rootPath);
else
    infoPath='';  % We assume that names include path
    rootPath=AddSlash(pwd);
end;
if ~iscell(miNames)
    miNames={miNames};
end;

%%
vmOk=0;
numErr=0;
maxErr=10; % Max number of error messages to print out

%% --------loop over files --------
for fileIndex=1:numel(miNames)
    %%
    disp(['Reading ' num2str(fileIndex) ' ' infoPath miNames{fileIndex}]);
    mi=ReadMiFile([infoPath miNames{fileIndex}]);
            originalAlpha=mi.ctf(1).alpha;
    if pars.modifiedAlpha>0
        mi.ctf(1).alpha=pars.modifiedAlpha;        
    end;
    if resetBasePath
        mi.basePath=rootPath;
    end;
    if ~isfield(mi,'mergeMode')
        mi.mergeMode=3;
    end;
    if ~(isfield(mi,'vesicle') && isfield(mi.vesicle,'x'))
        disp('No vesicles.');
        continue;  % skip on to the next file
    end;
    
    %     check if there is something to do.
    ampsNotRefined= size(mi.vesicle.s,2)<2 || all(all(mi.vesicle.s(:,2:end)==0));
    if numel(mi.vesicle.x)>0 && (ampsNotRefined || pars.overwrite) % something to do
        if ~isfield(mi,'log')
            mi.log=cell(0);
        end;
        
        %     Check to see if we have a generic vesicle model (central 5 points are
        %     essentially equal) in which case we try to load a good model.
        nvm=numel(mi.vesicleModel);
        vmCtr=ceil((nvm+1)/2);
        if (nvm < 3 || std(mi.vesicleModel(vmCtr-1:vmCtr+1))<1e-6...
                || forceNewModel)
            if ~vmOk  % Check only the first time through.
                mname=pars.modelMiName;
                pa='';
                if ~(exist(mname,'file') || batchMode)
                    disp('Getting an mi file for the vesicle membrane model');
                    [mname, pa]=uigetfile({'*mi.txt' '*mi.mat'},'Select an mi file for membrane model');
                end;
                if ischar(mname) % we have an mi file name
                    disp(['Loading the membrane model from ' pa mname]);                    vm=ReadMiFile([pa mname]);
                    vm=ReadMiFile(mname);
                    vmOk=1;
                    vmGood=1;
                end;
            end;
            if vmGood  % replace the model with the one in vm
                disp(['Using the membrane model from ' mname]);
                if vm.pixA == mi.pixA  % Check if the same pixel size
                    mi.vesicleModel=vm.vesicleModel;  % copy the model
                else                   % resample the model
                    disp('Resampling the model');
                    mi.vesicleModel=meDownsampleVesicleModel(...
                        vm.vesicleModel,mi.pixA/vm.pixA);
                end;
                %                 Replace the ampFactors too.
                for i=1:min(numel(vm.ctf),numel(mi.ctf))
                    mi.ctf(i).ampFactor=vm.ctf(i).ampFactor;
                end;
            else
                disp('Using the existing generic membrane model');
            end;
        else  % the existing model is a good one.
            disp('Using the existing  membrane model');
            mname=miNames{fileIndex};
            vmGood=1;
            vm=mi;
        end;
        
% ----------Read the image and normalize to fractional contrast----------
%         The final coordinates will all be with respect to the original
%         raw image, regardless of whether we are working from a merged
%         image or not. We mark this so:
        mi.useMicrographCoords=1;
%         Otherwise the coordinates would be, as before, with respect to
%         the padded "merged" image.

%           Load the full-sized, padded and scaled image, with M1
%           indicating the shift wrt the original micrograph.
        [m1,M1,ok,isRawImg]=meLoadNormalizedImage(mi,mi.padImageSize,'m');
        if ~ok
            disp(['No image found for ' mi.baseFilename]);
            numErr=numErr+1;
            if numErr>maxErr
                break;
            else
                continue; % Go on to the next file index.
            end;
        end;
%           Make downsampled copies for fitting. m4 is taken to be "full
%           size" for amplitude fitting. maxPixA ~ 3A.
        dsMin=pars.maxPixA/mi.pixA;
        ds4=NextNiceNumber(dsMin,5,1); % Downsampling for good fits
%         some possible values: 1, 2, 3, 4, 5, 6, 8, ....
%         ds4 is typically 2 or 3.
        disp(['First downsampling is by ' num2str(ds4)]);
        n4=round(mi.padImageSize/ds4);
        [m4,M4]=meDownsampleImage(m1,M1,n4);
        ds8=2*ds4;
%        The second downsampling is twice that. Images used for radius
%        fitting.
        n8=round(mi.padImageSize/ds8);
        [m8,M8]=meDownsampleImage(m4,M4,n8);

        % % %         pixAWork=pars.scl.ds0*mi.pixA;
        
%        % [m0, mergePath]=meReadMergedImage(mi);
%         [mergedName,ok]=CheckForAltImage([mi.procPath mi.baseFilename 'm.mrc'],sufExts);
%         if ok
%             disp(['Merged image: ' mergedName]);
%             m0=ReadEMFile(mergedName);
%         end;
%         if ok && numel(m0)>1 % we've got an image
%             n0=size(m0,1);  % merged image size
%             dsm=mi.imageSize(1)/size(m0,1);  % downsampling factor of merged image
%             if dsModel>dsm
%                 disp(['Merged image will be downsampled by ' num2str(dsModel/dsm) ' for fitting.']);
%                 m=DownsampleGeneral(m0,n0/(dsModel/dsm));
%             else
%                 m=m0;
%             end;
%             %     Check that we have a temp directory
%             if ~isfield(mi,'tempPath')
%                 mi.tempPath='Temp/';
%             end;
%             if ~exist(mi.tempPath,'dir')
%                 mkdir('Temp');
%             end;
%             if ~isfield(mi,'mask')
%                 mi.mask=[];
%             end;
%             
%             if pars.scaleOriginalAmplitudes~=1
%                 mi.vesicle.s=mi.vesicle.s*pars.scaleOriginalAmplitudes;
%             end;
%             
%             %                    mi.vesicle.s=mi.vesicle.s*dsm/2;  % scale down if not downsamping....
%             
            
            %%  -------------------All the work is done here---------------------------
            %   do nRounds fits; each time fit all the vesicles in the
            %   image.
            nRounds=numel(pars.fitModes);
            for ind=1:nRounds  % We loop through 2-3 times, first with no extra peaks,
%                                 and perhaps with no amplitude fitting
%                                 (except the constant term).
                %                 The first time through we fit no peaks because xPeakPositionA{1}=[]
                %                 later times include the extra peaks.
                p=struct;
%                 p.scl=pars.scl; % scaling of downsampled image
                p.extraPeaks=pars.xPeakPositionA{ind}/mi.pixA;
                p.extraSD=pars.xPeakSigmaA{ind}/mi.pixA;
                p.rTerms=pars.rTerms;
                p.disA=pars.disA;
                p.M8=M8;
                p.M4=M4;
                maxRTerms=numel(p.rTerms);
                % Set up parameters
                p.limitOrigNTerms=round(pars.fractionStartingTerms(ind)*maxRTerms+1);
                p.fitMode=pars.fitModes{ind};
                disp(['fitMode = ' p.fitMode]);
                p.radiusStepsA=pars.radiusStepsA;
                %                 Limit the number of amplitude terms.
                %                 On the first ind, use only 1-2 amp terms.  On later
                %                 iterations, use 2/3 times as many as radius terms.
                p.fracAmpTerms=pars.fractionAmpTerms(ind); % e.g. [.1 1].
                p.doPreSubtraction=pars.doPreSubtraction | ind>1; % force pre-subtraction in rounds
                p.listFits=pars.listFits;  % write out individual fit values?
                disp(['Round ' num2str(ind) '  ' miNames{fileIndex}]);
                %       --------- Do the fitting  -----------
                mi1=rsRefineVesicleFitsSub(mi,m8,m4,p,displayOn);
                mi=mi1;
            end;
            
            
            %% ---------------Outputting------------------
%             Write the mi file.
            outName='';
            mi.ctf(1).alpha=originalAlpha;
            if writeMiFile
                mi.log{end+1,1}=['rsRefineVesicleFits ' TimeStamp];
                outName=WriteMiFile(mi,[infoPath miNames{fileIndex}]);
                disp([infoPath outName ' saved']);
            end;
            
            %             Compute and store model vesicles
            pixAModel=mi.pixA*dsModel;  % pixel size of working merged image
%             pixA0=mi.pixA*dsm;
            
            dfc=.1;  % display filter relative to raw data
            %%
            if displayOn
                figure(1); clf;
                imags(GaussFilt(m1,dfc*dsModel));  % show the unsubtracted image
                title(['Original image ' outName],'interpreter','none');
                drawnow;
            end;
            
%             We'll make the vesicles at the m4 size and scale up.
            scl4=struct;
            scl4.n=size(m4);
            scl4.M=M4;
            disp('Making the final vesicle models');
            vs4=meMakeModelVesicles(mi,scl4,find(mi.vesicle.ok(:,3)));
            vs1=Crop(Downsample(vs4,M4(1,1)*size(vs4)),size(m1));  % scale up if needed to match m0
            
            if displayOn
                imags(GaussFilt(m4-vs4,dfc*M4(1,1)));
                title(['Final subtraction ' outName],'interpreter','none');
                drawnow;
            end;
            
%             
%             if pars.writeVesFiles  % Write .mrc and .jpg files.
%                 outVesName=[mi.tempPath mi.baseFilename 'v'];
%                 %             WriteMRC(vsm,pixA0,[outVesName '.mrc']);
%                 %             WriteJpeg(vsm,outVesName);
%                 WriteMRC(vs4,pixAModel,[outVesName '.mrc']);
%                 WriteJpeg(vs4,outVesName,0);
%                 %             imwrite(uint8(imscale(rot90(vs1),256,0)),[outVesName '.jpg']);
%                 disp([outVesName ' saved']);
%             end;
%             
            if ~isfield(mi,'procPath_sm')
                procPath_sm=mi.procPath;
            else
                procPath_sm=mi.procPath_sm;
            end;
            
            if pars.writeSubMRC  % write an MRC file
                if isRawImg % our input is a raw micrograph, make the output the same size.
                    mSub=Crop(m1-vs1,mi.imageSize);
                    outSubName=[mi.procPath mi.baseFilename '_v.mrc'];
                else
                    mSub=m1-vs1;
                    outSubName=[mi.procPath mi.baseFilename 'mv.mrc'];
                end;
                WriteMRC(mSub,mi.pixA,outSubName);
                disp([outSubName ' saved']);
            end;

            smSize=round(size(m1)/pars.dsSmall);
            if pars.writeSmallMRC
                outSmallName=[procPath_sm mi.baseFilename 'ms.mrc'];
                ms=Downsample(m1,smSize);
                WriteMRC(ms,mi.pixA*pars.dsSmall,outSmallName);
                disp([outName ' saved.']);
            end;
            if pars.writeSmallSubMRC
                outSubName=[procPath_sm mi.baseFilename 'mvs.mrc'];
                mvs=Downsample(m1-vs1,smSize);
                WriteMRC(mvs,mi.pixA*pars.dsSmall,outSubName);
                disp([outSubName ' saved, ' num2str(smSize) ' pixels']);
            end;
            
    else  % No vesicles have been found to refine
        if numel(mi.vesicle.x)<1
            disp('  ...no vesicles');
        else
            disp('  ...already refined.');
        end;
    end;
    disp(' ');
end; % for fileIndex
disp('Done.');
disp(' ');

