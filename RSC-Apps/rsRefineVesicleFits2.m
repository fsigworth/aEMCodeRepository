function rsRefineVesicleFits2(miNames,mpars)
% function rsRefineVesicleFits2(miNames,mpars)
% Find and subtract vesicles
% and return the updated mi containing the vesicle coordinates, radius
% expansion and amplitude angular expansion.
%
% Cleaned-up version. Fits distorted vesicles.
% miNames is a cell array of mi file names (with path)
% mpars is a struct of parameters (see the default values set below)
% If desired, a vesicle model file is stored in Temp/ as <basename>v.mrc.
% If desired, a subtracted file is stored in Merged/ as <basename>mv.mrc.

% if nargin<1
%     miNames=[];
% end;
disp('rsRefineVesicleFits2 startup.');


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
dpars.writeJpeg=1;
dpars.writeSmallSubMRC=1; % Write out a downsampled subtracted image
dpars.writeSubJpeg=1;
dpars.jpegPath='Merged_jpeg/';
dpars.dsSmall=4; % downsampling factor for small output images
dpars.maxPixA=6.5;  % downsampled image resolution for radius fitting
dpars.forceDs4=4;  % Or, if nonzero, use this fixed downsampling factor instead
% dpars.writeVesFiles=0;   % Write vesicle models into Temp/

dpars.resetBasePath=1;   % update the basePath field of the mi file to the current directory.
dpars.modifiedAlpha=.02;  % Good alpha value for lipids.
dpars.modifiedB=.1;    % nonzero to change mi.CTF.B
dpars.modifiedImageNormScale=.2;
dpars.maxVesiclesToFit=inf;
% Number of terms (for both radius and amplitude fitting) is set thusly:
%    nTerms=find(vesicle.r(ind,1) < pars.rTerms/mi.pixA,1);
% i.e. nTerms is the index of last entry of rTerms smaller than our radius.
% dpars.rTerms=[100 150 200 250 250 300 inf];
dpars.rTerms=[100 150 200 300 inf];  % quick
% dpars.rTerms=[150 200 inf];  % for Mengqiu
dpars.limitOrigNTerms=4; % max number of radius terms to fit first
dpars.stepNTerms=3;
% Define the fitting mode for each round
% dpars.fitModes={'RadiusOnly' 'RadiusAndLin'};
dpars.fitModes={'RadiusOnly' 'LinOnly'}; % First fit shape, then amplitudes
% dpars.fitModes={'LinOnly'};
%    dpars.fitModes={'RadiusAndLin'};
dpars.fractionAmpTerms=[0 1]; % fraction of amp terms to use
%  To avoid crazy fits, we repeat the whole radius-only fitting with the base
%  radius perturbed by these steps (in angstroms) and pick the best.
dpars.radiusStepsA=[-100 -50 0 50]; % repeat radius-only fitting with perturbed r(1)
dpars.disA=1200;  % size of the fit window, in angstroms
% dpars.disA=1600;  % for Mengqiu

%  We add extra peaks in the scattering profile when fitting amplitudes
dpars.peakPositionA=[-37 0 37];  % empirical default.  Works a bit better than [37 37]
dpars.xPeakSigmaA={5 5}; % width of extra Gaussian peaks, in angstrom
% There must be as many elements as rounds.
dpars.suffix='_v2+peaks+300iters_conv0.1_newModel_';

% -----------Merge the defaults with the given mpars-----------
pars=SetOptionValues(dpars,mpars,1);

%     The following must have at least as many elements as dpars.fitModes!
pars.xPeakPositionA=cell(1,numel(pars.fitModes));
pars.xPeakPositionA{end}=pars.peakPositionA; % add the peaks to the last round.

% -----------------------

% See if we are on the cluster.
host=getenv('HOSTNAME');
disp(['host name: ' host]);
isCluster=any(strncmp(host,{'c' 'p'},1)); % if we are on a farnam node
batchMode=isCluster  && (exist('miNames','var') && numel(miNames)>0);

displayOn=~batchMode;
% displayOn=0; % 1: show results at end of each fit; 2 update while fitting.

% dsModel=2;  % net downsampling of "full res" model relative to orig micrograph

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
    tic; % We'll time the processing of this micrograph
    disp(['Reading ' num2str(fileIndex) ' ' infoPath miNames{fileIndex}]);
    mi=ReadMiFile([infoPath miNames{fileIndex}]);

    % Make changes to the mi structure
    originalAlpha=mi.ctf(1).alpha;
    if pars.modifiedAlpha>0
        mi.ctf(1).alpha=pars.modifiedAlpha;        
    end;
    if pars.modifiedB>0
        mi.ctf(1).B=pars.modifiedB;
    end;
    if pars.resetBasePath
        mi.basePath=rootPath;
    end;
    if ~isfield(mi,'mergeMode')
        mi.mergeMode=3;
    end;

    % No vesicle entries?
    if ~(isfield(mi,'vesicle') && isfield(mi.vesicle,'x'))
        disp('No vesicles.');
        continue;  % skip on to the next file
    end;
    
    %     Check if there is something to do.
    ampsNotRefined= size(mi.vesicle.s,2)<2 || all(all(mi.vesicle.s(:,2:end)==0));
    if numel(mi.vesicle.x)>0 && (ampsNotRefined || pars.overwrite) % something to do
        if ~isfield(mi,'log')
            mi.log=cell(0);
        end;
        
        %     Check to see if we have a generic vesicle model (central 5 points are
        %     essentially equal) in which case we try to load a good model.
        nvm=numel(mi.vesicleModel);
        vmCtr=ceil((nvm+1)/2);
        if (nvm < 3 || std(mi.vesicleModel(vmCtr-1:vmCtr+1))<1e-6)
            disp('--only has a generic vesicle model.')
        end;

% ----------Read the image and normalize to fractional contrast----------
%         The final coordinates will all be with respect to the original
%         raw image, regardless of whether we are working from a merged
%         image or not. We mark this so:
        mi.useMicrographCoords=1;
%         Otherwise the coordinates would be, as before, with respect to
%         the padded "merged" image.

%           Load the full-sized, padded and scaled image, with matrix M1
%           indicating the shift wrt the original micrograph.
        [m1,M1,ok,isRawImg]=meLoadNormalizedImage(mi,mi.padImageSize,'m'); 
%         ---We assume m1 has the pixel size mi.pixA
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
if pars.forceDs4
    ds4=pars.forceDs4;
else
        dsMin=pars.maxPixA/mi.pixA;
        ds4=NextNiceNumber(dsMin,5,1); % Downsampling for good fits
%         some possible values: 1, 2, 3, 4, 5, 6, 8, ....
%         ds4 is typically 2 or 3.
end;
        disp(['First downsampling is by ' num2str(ds4)]);
        n4=round(mi.padImageSize/ds4);
        [m4,M4]=meDownsampleImage(m1,M1,n4);
%        The second downsampling is twice that. This is used for radius
%        fitting.
        ds8=2*ds4;
        n8=round(mi.padImageSize/ds8);
        [m8,M8]=meDownsampleImage(m4,M4,n8);

            
            %%  -------------------All the work is done here---------------------------

            nRounds=numel(pars.fitModes);
            for iRound=1:nRounds  % We loop through 2-3 times, first with no extra peaks,
%                                 and perhaps with no amplitude fitting
%                                 (except the constant term).
                %                 The first time through we fit no peaks because xPeakPositionA{1}=[]
                %                 later times include the extra peaks.
                p=struct;
                p.extraPeaks=pars.xPeakPositionA{iRound}/mi.pixA; % will be [] if no peaks applied.
                p.extraSD=pars.xPeakSigmaA{iRound}/mi.pixA;
                p.rTerms=pars.rTerms;
                p.maxVesiclesToFit=pars.maxVesiclesToFit;
%                 if p.maxVesiclesToFit< size(mi.vesicle.ok,1)  %%%%%%%kill the other vesicles....
%                     mi.vesicle.ok(p.maxVesiclesToFit:end,1)=0;
%                 end;
                p.disA=pars.disA; % fit window size, in A.
                p.M8=M8;
                p.M4=M4;
%                 maxRTerms=numel(p.rTerms);
                % Set up parameters
                p.limitOrigNTerms=pars.limitOrigNTerms;
                p.fitMode=pars.fitModes{iRound};
                disp(['fitMode = ' p.fitMode]);
                p.radiusStepsA=pars.radiusStepsA; % perturbation of radii in separate fits.
                %                 Limit the number of amplitude terms.
                %                 On the first ind, use only 1-2 amp terms.  On later
                %                 iterations, use 2/3 times as many as radius terms.
                p.fracAmpTerms=pars.fractionAmpTerms(iRound); % e.g. [0 1].
                p.doPreSubtraction=pars.doPreSubtraction | iRound>1; % force pre-subtraction in rounds
                p.listFits=pars.listFits;  % write out individual fit values?
                disp(['Round ' num2str(iRound) '  ' miNames{fileIndex}]);
                %       --------- Do the fitting  -----------
                mi1=rsRefineVesicleFitsSub(mi,m8,m4,p,displayOn);
                mi=mi1;
            end;
            
            
            %% ---------------Outputting------------------

            dfc=.1;  % display filter in A^-1
            %%
            outName=miNames{fileIndex};
            if displayOn
                figure(1); clf;
                imags(GaussFilt(m1,dfc*mi.pixA));  % show the unsubtracted image
                title(['Original image ' outName],'interpreter','none');
                drawnow;
            end;
            
%             Compute and store model vesicles
%             We'll make the vesicles at the m4 size and scale up.
            scl4=struct;
            scl4.n=size(m4);
            scl4.M=M4; % we make it the same size as 1/ds4, typically 1/4
            disp('Making the final vesicle models');
            vs4=meMakeModelVesicles(mi,scl4,find(mi.vesicle.ok(:,3)));
            vs1=Crop(Downsample(vs4,M4(1,1)*size(vs4)),size(m1));  % scale up if needed to match m0
            
            if displayOn
                imags(GaussFilt(m4-vs4,dfc*M4(1,1)));
                title(['Final subtraction ' outName],'interpreter','none');
                drawnow;
            end;
            
%             %             Write the mi file.
            mi.ctf(1).alpha=originalAlpha;
            if pars.writeMiFile
                mi.log{end+1,1}=['rsRefineVesicleFits ' TimeStamp];
                outName=WriteMiFile(mi,[infoPath miNames{fileIndex}]);
                disp([infoPath outName ' saved']);
            end;
 
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
            
            suffix=pars.suffix

            if pars.writeSubMRC  % write an MRC file
                CheckAndMakeDir(mi.procPath);
                if isRawImg % our input is a raw micrograph, make the output the same size.
                    % ...and undo the image normalization so it matches the raw
                    % image.
                    mSub=(Crop(m1-vs1,mi.imageSize)/mi.imageNormScale)+mi.imageMedian;
                    outSubName=[mi.procPath mi.baseFilename suffix '_v.mrc'];
                else
                    mSub=m1-vs1;
                    outSubName=[mi.procPath mi.baseFilename suffix 'mv.mrc'];
                end;
                WriteMRC(mSub,mi.pixA,outSubName);
                disp([outSubName ' saved']);
            end;

            if pars.writeSmallMRC || pars.writeJpeg
            smSize=round(size(m1)/pars.dsSmall);
                ms=Downsample(m1,smSize);
            end;
            
            if pars.writeSmallMRC
                CheckAndMakeDir(procPath_sm);
                outSmallName=[procPath_sm mi.baseFilename suffix 'ms.mrc'];
                WriteMRC(ms,mi.pixA*pars.dsSmall,outSmallName);
                disp([outName ' saved.']);
            end;
            if pars.writeJpeg
                outJpegName=[pars.jpegPath mi.baseFilename suffix 'ms.jpg'];
                CheckAndMakeDir(pars.jpegPath);
                WriteJpeg(ms,outJpegName);
                disp([outJpegName ' saved.']);
            end;
            
            if pars.writeSmallSubMRC || pars.writeSubJpeg
                mvs=Downsample(m1-vs1,smSize);
            end;
            if pars.writeSubJpeg
                CheckAndMakeDir(pars.jpegPath);
                outJpegName=[pars.jpegPath mi.baseFilename suffix 'mvs.jpg'];
                WriteJpeg(mvs,outJpegName);
                disp([outJpegName ' saved.']);
            end;
        toc; % display the time used for this refined micrograph.
            
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

