function rsRefineVesicleFits(miNames,mpars)
% function rsRefineVesicleFits(miNames,mpars)
% Fits distorted vesicles.
% Revised version with acceleration options.
% Given an info structure mi, find and subtract vesicles
% and return the updated mi containing the vesicle coordinates, radius
% expansion and amplitude angular expansion.
% If desired, a vesicle model file is stored in Temp/ as <basename>v.mrc.
% If desired, a subtracted file is stored in Merged/ as <basename>mvz.tif.

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
dpars.doPreSubtraction=1;
dpars.listFits=1;  % print out each fit's parameters
dpars.scaleOriginalAmplitudes=1;
dpars.scaleOriginalAmplitudes=0.5; %%%%%%
% Number of terms (for both radius and amplitude fitting) is set thusly:
%    nTerms=find(vesicle.r(ind,1) < pars.rTerms/mi.pixA,1);
% i.e. nTerms is the index of last entry of rTerms smaller than our radius.
dpars.rTerms=[100 150 200 250 250 300 inf];

% Define the fitting mode for each round
% dpars.fitModes={'RadiusOnly' 'RadiusAndLin'};
    dpars.fitModes={'RadiusOnly' 'LinOnly'};
%    dpars.fitModes={'RadiusAndLin'};
dpars.fractionStartingTerms=[.5 1]; % total terms to use in each round
dpars.fractionAmpTerms=[0 1];
dpars.radiusStepsA=[-100 -50 0 50]; % repeat radius-only fitting with perturbed r(1)
% Extra peaks in the scattering profile
dpars.peakPositionA=[-37 0 37];  % empirical default.  Works a bit better than [37 37]
dpars.targetPixA=10;  % downsampled image resolution for radius fitting

dpars.xPeakSigmaA={5 5}; % width of extra Gaussian peaks, in angstrom

% Merge the defaults with the given mpars
pars=SetOptionValues(dpars,mpars,1);

% pars

%     The following must have at least as many elements as dpars.fitModes!
pars.xPeakPositionA=cell(1,numel(dpars.fitModes));
pars.xPeakPositionA{end}=pars.peakPositionA;
% pars.xPeakPositionA={[] pars.peakPositionA}; % centers of extra peaks, in angstrom
% -----------------------

% See if we are on the cluster.
host=getenv('HOSTNAME');
disp(['host name: ' host]);
isCluster=strncmp(host,'c',1); % if we are on a farnam node
batchMode=isCluster  && (exist('miNames','var') && numel(miNames)>0);

displayOn=~batchMode;
% displayOn=0; %####### 1: show results at end of each fit; 2 update while fitting.

dsModel=4;  % net downsampling of "full res" model relative to orig micrograph
dsSmall=4;  % net downsampling of 'small' output file.

forceNewModel=0;   % Always ask the user to select a new refined model
  % on the first micrograph (can be from the same micrograph)
resetBasePath=1;   % update the basePath field of the mi file.

writeMiFile=pars.writeMiFile;     % Save the updated mi file

writeSubZTiff=0;    % Write subtracted images into Merged/
suffix='';  % extra character ? for XXmv?z.tif output file
writeSubMRC=1;
writeSmallSubMRC=1;
writeVesFiles=0;   % Write vesicle models into Temp/

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
    infoPath='Info/'; %%%%%%%%%%%%%%
    infoPath='';
    rootPath=AddSlash(pwd);
end;
if ~iscell(miNames)
    miNames={miNames};
end;

%%
vmOk=0;

%% --------loop over files --------
for fileIndex=1:numel(miNames)
    %%
    disp(['Reading ' infoPath miNames{fileIndex}]);
    mi=ReadMiFile([infoPath miNames{fileIndex}]);
    if resetBasePath
        mi.basePath=rootPath;
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
        %     essentially equal)
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
        
        % Read the image and normalize to fractional contrast
        sufExts={'.mrc' 's.mrc' 'z.tif'}; % suffix and extension options for small, compressed or full files.
        % drawnow;
        % [m0, mergePath]=meReadMergedImage(mi);
        [mergedName,ok]=CheckForAltImage([mi.procPath mi.baseFilename 'm.mrc'],sufExts);
        if ok
            disp(['Merged image: ' mergedName]);
            m0=ReadEMFile(mergedName);
        end;
        if ok && numel(m0)>1 % we've got an image
            n0=size(m0,1);  % merged image size
            dsm=mi.imageSize(1)/size(m0,1);  % downsampling factor of merged image
            if dsModel>dsm
                disp(['Merged image will be downsampled by ' num2str(dsModel/dsm) ' for fitting.']);
                m=DownsampleGeneral(m0,n0/(dsModel/dsm));
            else
                m=m0;
            end;
            %     Check that we have a temp directory
            if ~isfield(mi,'tempPath')
                mi.tempPath='Temp/';
            end;
            if ~exist(mi.tempPath,'dir')
                mkdir('Temp');
            end;
            if ~isfield(mi,'mask')
                mi.mask=[];
            end;
            
            if pars.scaleOriginalAmplitudes~=1
                mi.vesicle.s=mi.vesicle.s*pars.scaleOriginalAmplitudes;
            end;
            
            %                    mi.vesicle.s=mi.vesicle.s*dsm/2;  % scale down if not downsamping....
            
            
            %%  -------------------All the work is done here---------------------------
            %   do nRounds fits; each time fit all the vesicles in the
            %   image.
            nRounds=numel(pars.fitModes);
            for ind=1:nRounds  % We loop through 2-3 times, first with no extra peaks
                %                 The first time through we fit no peaks because xPeakPositionA{1}=[]
                %                 later times include the extra peaks.
                p=struct;
                p.extraPeaks=pars.xPeakPositionA{ind}/mi.pixA;
                p.extraSD=pars.xPeakSigmaA{ind}/mi.pixA;
                p.rTerms=pars.rTerms;
                p.targetPixA=pars.targetPixA;
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
                %                   Fit all the vesicles in this mi file.
                mi1=rsRefineVesicleFitsSub(mi,m,p,displayOn);
                mi=mi1;
            end;
            
            %% Outputting
            outName='';
            
            if writeMiFile
                mi.log{end+1,1}=['rsRefineVesicleFits ' TimeStamp];
                outName=WriteMiFile(mi,miNames{fileIndex});
                disp([outName ' saved']);
            end;
            
            %             Compute and store model vesicles
            pixAModel=mi.pixA*dsModel;  % pixel size of working merged image
            pixA0=mi.pixA*dsm;
            
            dfc=.1;  % display filter relative to raw data
            %%
            if displayOn
                figure(1); clf;
                imags(GaussFilt(m,dfc*dsModel));  % show the unsubtracted image
                title(['Original image ' outName],'interpreter','none');
                drawnow;
            end;
            
            disp('Making the final vesicle models');
            vs1=meMakeModelVesicles(mi,size(m),find(mi.vesicle.ok(:,3)));
            vsm=Downsample(vs1,n0);  % scale up if needed to match m0
            
            if displayOn
                imags(GaussFilt(m0-vsm,dfc*dsm));
                title(['Final subtraction ' outName],'interpreter','none');
                drawnow;
            end;
            
            
            if writeVesFiles  % Write .mrc and .jpg files.
                outVesName=[mi.tempPath mi.baseFilename 'v'];
                %             WriteMRC(vsm,pixA0,[outVesName '.mrc']);
                %             WriteJpeg(vsm,outVesName);
                WriteMRC(vs1,pixAModel,[outVesName '.mrc']);
                WriteJpeg(vs1,outVesName,0);
                %             imwrite(uint8(imscale(rot90(vs1),256,0)),[outVesName '.jpg']);
                disp([outVesName ' saved']);
            end;
            
            if writeSubZTiff  % write a zTiff file
                ztPars=struct;
                ztPars.lfCutoff=.1;
                ztPars.snrRatio=300;
                outSubName=[mi.procPath mi.baseFilename 'mv' suffix 'z.tif'];
                WriteZTiff(m0-vsm,pixA0,outSubName,ztPars);
                disp([outSubName ' saved']);
            end;
            
            if writeSubMRC  % write an MRC file
                outSubName=[mi.procPath mi.baseFilename 'mv' suffix '.mrc'];
                WriteMRC(m0-vsm,pixA0,outSubName);
                disp([outSubName ' saved']);
            end;
            if writeSmallSubMRC
                nSmall=mi.imageSize/dsSmall;
                outSubName=[mi.procPath mi.baseFilename 'mvs' suffix '.mrc'];
                WriteMRC(Downsample(m0-vsm,nSmall),mi.pixA*dsSmall,outSubName);
                disp([outSubName ' saved']);
            end;
        else
            disp('  ...no image read.');
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

