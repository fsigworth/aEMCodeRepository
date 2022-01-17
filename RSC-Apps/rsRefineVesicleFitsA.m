function rsRefineVesicleFitsA(miNames,mpars)
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
% pa=AddSlash(fileparts(which('rsRefineVesicleFits'))); % Get our directory
dpars=struct;
dpars.overwrite=1; % Fit even if a fit was done before.
dpars.skipFitting=0; % don't fit at all, just make and write image files
dpars.writeMiFile=1;
dpars.doPreSubtraction=1; % subtract all vesicles except the one being fitted
dpars.listFits=1;  % print out each fit's parameters
dpars.maxVesiclesToFit=inf;
dpars.writeSmallMRC=1; % Write out image downsampled to M4. Otherwise, we
dpars.writeJpeg=1;
dpars.writeSubMRC=1; % Write out subtracted image <basename>mv.mrc
dpars.writeSmallSubMRC=1; % Write out a downsampled subtracted image
dpars.writeSubJpeg=1;
dpars.jpegPath='Merged_jpeg/';

dpars.dsSmall=4; % downsampling factor for small output images
dpars.dsSmaller=8;

dpars.fHP=.003; % Highpass filter for fitting
dpars.dfc=.05;  % display filter in A^-1

% dpars.writeVesFiles=0;   % Write vesicle models into Temp/

dpars.modifiedAlpha=.02;  % Good alpha value for lipids.
dpars.modifiedB=0.1;    % nonzero to change mi.CTF.B
% Number of terms (for both radius and amplitude fitting) is set thusly:
%    nTerms=find(vesicle.r(ind,1) < pars.rTerms/mi.pixA,1);
%    i.e. nTerms is the index of last entry of rTerms smaller than our radius.
%    Note that 2 terms is the same as 1 term.
%    e.g. dpars.rTerms=[100 100 120 150 200 300 inf];
%    in this example, r>100A gets 3 terms, r>120 gets 4 terms, etc.
% dpars.rTerms=[100 100 140 180 220 250 300 inf]; % set 1
dpars.rTerms=[100 100 150 200 250 300 inf]; % set 2 good for Kv
% dpars.rTerms=[120 120 180 240 300 inf]; % set 3
dpars.aTerms=dpars.rTerms;

dpars.stepNTerms=2;
dpars.preservedTerms=3; % how many radius terms to preserve from previous fitting
dpars.maxMaskLayers=inf;

% Define the fitting mode for each round
dpars.fitRadiusModes=[1 0]; % Only one pass, fit radius (shape) alone.
% dpars.fitRadiusModes=[1 0]; % First fit shape, then amplitudes

%  To avoid crazy fits, we repeat the whole radius-only fitting with the base
%  radius perturbed by these steps (in angstroms) and pick the best.
% dpars.radiusStepsA=[-100 -50 0 50]; % repeat radius-only fitting with perturbed r(1)
dpars.radiusStepsA=[-30]; % this value cancels the initial step magnitude.

dpars.disA=1200;  % size of the fit window, in angstroms
dpars.nIters=300; % basic number of Simplex iterations
dpars.fHP=.003; % highpass for fitting, in A^-1
dpars.displayOn=2;
%  We can add extra peaks in the scattering profile when fitting amplitudes
dpars.peakPositionA=[]; % for Krios Kv data, no peak! Old favorite [-37 0 37]
dpars.peakSigmaA=5; % width of extra Gaussian peaks, in angstrom
dpars.extraRound=0;
% -----------Merge the defaults with the given mpars-----------
pars=SetOptionValues(dpars,mpars);

% -----------------------

% % See if we are on the cluster.
% host=getenv('HOSTNAME');
% disp(['host name: ' host]);
% isCluster=any(strncmp(host,{'c' 'n'},1)); % 1st char of name if we are on a farnam node

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
end;
if ~iscell(miNames)
    miNames={miNames};
end;

%%
numErr=0;
maxErr=10; % Max number of error messages to print out

%% --------loop over files --------
for fileIndex=1:numel(miNames)
    %%
    disp(['Reading ' num2str(fileIndex) ' ' infoPath miNames{fileIndex}]);
    mi=ReadMiFile([infoPath miNames{fileIndex}]);
    % Make changes to the mi structure
    originalAlpha=mi.ctf(1).alpha;
    if pars.modifiedAlpha>0
        mi.ctf(1).alpha=pars.modifiedAlpha;
    end;
    originalB=mi.ctf(1).B;
    if pars.modifiedB>0
        mi.ctf(1).B=pars.modifiedB;
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
        ds4=pars.dsSmall;
        disp(['First downsampling is by ' num2str(ds4)]);
        n4=round(mi.padImageSize/ds4);
        [m4,M4]=meDownsampleImage(m1,M1,n4);
        pars.M4=M4;
        ds8=pars.dsSmaller;
        %        The second downsampling is twice that. Images used for radius
        %        fitting.
        n8=round(mi.padImageSize/ds8);
        [m8,M8]=meDownsampleImage(m4,M4,n8);
        pars.M8=M8;
        if ~pars.skipFitting

            %%  -------------------All the work is done here---------------------------
            %   do nRounds fits; each time fit all the vesicles in the
            %   image.
            figure(2); % this is where we show the progress
            nRounds=numel(pars.fitRadiusModes);
            for ind=1:nRounds  % We loop through 2-3 times, first with no extra peaks,
                disp(['Round ' num2str(ind) '  ' miNames{fileIndex}]);
                %                                 and perhaps with no amplitude fitting
                %                                 (except the constant term).
                %                 The first time through we fit no peaks because xPeakPositionA{1}=[]
                %                 later times include the extra peaks.
                p=pars;
                %                 p.scl=pars.scl; % scaling of downsampled image
                p.rTerms=pars.rTerms;
                % Set up parameters
                %                 p.limitOrigNTerms=round(pars.fractionStartingTerms(ind)*maxRTerms+1);
                p.fitRadius=pars.fitRadiusModes(ind);
                disp(['fitRadius = ' num2str(p.fitRadius)]);
                if p.fitRadius
                    p.M=pars.M8; % smaller size for radius
                    mWork=m8;
                    p.radiusStepsA=pars.radiusStepsA;
                else % fitting Amplitudes
                    p.M=pars.M4;
                    mWork=m4;
                end;
                %       --------- Do the fitting  -----------
                mi1=rsRefineVesicleFitsSubA(mi,mWork,p);
                mi=mi1;
            end; % for rounds

            %% --------------Write the mi file.-----------------
            %
            outName='';
            if pars.writeMiFile
                mi.ctf(1).alpha=originalAlpha; % restore the original values
                mi.ctf(1).B=originalB;
                mi.log{end+1,1}=['rsRefineVesicleFits ' TimeStamp];
                outName=WriteMiFile(mi,[infoPath miNames{fileIndex}]);
                disp([infoPath outName ' saved']);
            end;
        else
            mi1=mi;
        end; % ~skipFitting


        %% ----------display the result -------
        if pars.displayOn
            figure(1); clf;
            imags(GaussFilt(m1,pars.dfc*mi.pixA));  % show the unsubtracted image
            title(['Original image ' mi.baseFilename],'interpreter','none');
            drawnow;
        end;

        %             We'll make the vesicles at the m4 size and scale up.
        scl4=struct;
        scl4.n=size(m4);
        scl4.M=M4;
        ds4=M4(1,1);

        disp('Making the final vesicle models');
        vs4=meMakeModelVesicles(mi1,scl4,find(mi.vesicle.ok(:,3)),1,0);
        % Use mi1 for the modified CTF.
        %             Upsample to match m1
        vs1=Crop(Downsample(vs4,M4(1,1)*size(vs4)),size(m1));  % scale up and shift to match m1

        if pars.displayOn
            figure(3); clf;
            imags(GaussFilt(m1-vs1,pars.dfc*mi.pixA));
            title(['Final subtraction ' mi.baseFilename],'interpreter','none');
            drawnow;
        end;

        %% -------------write out the images-----------
        if ~isfield(mi,'procPath_sm')
            % construct it!
            [pa,nm]=fileparts(mi.procPath); % Trick to get rid of traling filesep
            mi.procPath_sm=[pa nm '_sm/'];
        end;

        if pars.writeSubMRC  % write an MRC file
            CheckAndMakeDir(mi.procPath);
            if isRawImg % our input is a raw micrograph, make the output the same size.
                % ...and undo the image normalization so it matches the raw
                % image.
                mSub=(Crop(m1-vs1,mi.imageSize)/mi.imageNormScale)+mi.imageMedian;
                outSubName=[mi.procPath mi.baseFilename '_v.mrc'];
            else
                mSub=m1-vs1;
                outSubName=[mi.procPath mi.baseFilename 'mv.mrc'];
            end;
            WriteMRC(mSub,mi.pixA,outSubName);
            disp([outSubName ' saved']);
        end;

        if pars.writeSmallMRC || pars.writeJpeg
            smSize=mi.padImageSize/ds4;
            ms=Downsample(m1,smSize);
        end;

        if pars.writeSmallMRC
            CheckAndMakeDir(mi.procPath_sm);
            outSmallName=[mi.procPath_sm mi.baseFilename 'ms.mrc'];
            WriteMRC(ms,mi.pixA*ds4,outSmallName);
            disp([outSmallName ' saved.']);
        end;
        if pars.writeJpeg
            outJpegName=[pars.jpegPath mi.baseFilename 'ms.jpg'];
            CheckAndMakeDir(pars.jpegPath);
            WriteJpeg(ms,outJpegName);
            disp([outJpegName ' saved.']);
        end;

        if pars.writeSmallSubMRC || pars.writeSubJpeg
            smSize=mi.padImageSize/ds4;
            mvs=Downsample(m1-vs1,smSize);
        end;

        if pars.writeSmallSubMRC
            ds8=pars.M8(1,1);
            CheckAndMakeDir(mi.procPath_sm);
            outSmallSubMRCName=[mi.procPath_sm mi.baseFilename 'mvs.mrc'];
            WriteMRC(mvs,mi.pixA*ds8,outSmallSubMRCName);
            disp([outSmallSubMRCName ' saved.']);
        end;
        if pars.writeSubJpeg
            CheckAndMakeDir(pars.jpegPath);
            outJpegName=[pars.jpegPath mi.baseFilename 'mvs.jpg'];
            WriteJpeg(mvs,outJpegName);
            disp([outJpegName ' saved.']);
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

