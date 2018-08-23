function [ri,si,moi]=reMakeRunInfo
% Sets up the structures ri, si and moi.
% ri is the run info for all workers; si is the entire raw stack; moi is
% the initial model. We assume that we in a directory containing
% Reconstructions/ and Stack/ .  Within Reconstructions/ the main program
% will create directories such as 96a/ where 'a' is given by pathSuffix.

mode2D=1;

% isSmall=false;  % run a small problem, e.g. on laptop
% isSmall=true;
isSmall=0;
isShort=1;
normalNTwins=2;  % for 3D only.

pathSuffix='s10v1r1';

% ---- Flags ----
fl=struct;
fl.mode2D=mode2D;
fl.useRandomRefs=2;  % noise only
fl.makeFakeData=false;
fl.doShowPriors=false;
fl.makeRawProbs=false;
fl.allowRestart=true;
fl.restartRoiCollection=true; %%%%
fl.useAltImgs=true;  % use the unsubtracted stack for alternative reconstruction.
fl.useAltImgs=~isShort; %%% dont use them for short problem.
fl.removeRings=false;


fl.writeAllGroupFiles=false;
fl.saveRoi=true;  % always save the composite roi structure for each twin
fl.saveMoi=true;   % save the composite moi structures
fl.saveFigures=true;
fl.useParFor=false;

% ---------Create the initial model---------
moi=struct;  % initialize the intensive model variables
moi.sigmaN=1;  % start with a high value
% moi.sigmaN=.1;  %%%%%%%%% start with a high value
moi.sigmaC=1;
moi.sigmaG=1;
moi.imgAmp=.001;
moi.b0=0;
moi.pRefs=1;
moi.refVols=0;  % placeholder
moi.pVols=1;    % placeholder
moi.intensiveFields=fieldnames(moi);

% ----------Create the run info----------
ri=struct;
ri.flags=fl;
ri.usePrior=1;
%ri.startMapName='KvMap.mat';
ri.startMapName='KvMapSkew.mat';
ri.startRes=30;

ri.thresholds=[1e-4 1e-4 1e-4]*.01;  % Thresh for trans, alpha, refs
%ri.thresholds=[1e-4 1e-4 1e-4]*0;  % Thresh for trans, alpha, refs
ri.refScale=.1;     % reference amplitude factor
ri.symmetry=4;
%ri.volMaskType='sphere';
ri.volMaskType='sphere+ellipsoid';
%ri.volMaskType='cylinder';

% Number of volumes
ri.nVols=1;
moi.pVols=ones(ri.nVols,1);

% ------Angle and shift search parameters, and start assigning the run info
% ri----------
ri.isos=[0 1]; % to handle rso and iso particles.
ri.maxShiftU=6; % Maximum translation searched
%ri.maxShiftU=4; % Maximum translation searched %%%%%%

ri.maxTxFraction=.8;

ri.angleStepsU=[2 2 2];  % alpha, beta, gamma steps in degrees
ri.angleStepsU=[2 10 10];  % alpha, beta, gamma steps in degrees

% if isSmall || isShort
%    ri.maxShift=4;
%    ri.angleSteps=[6 6 10];
% end;

ri.angleLimits=[20 0 360/ri.symmetry];  % +-alpha range;  beta start; max gamma.

ri.startIter=1;
ri.nIters=30;

% Final image sizes
ri.nCropU=120;  % We first crop the input images to this size
ri.nFinal=ri.nCropU;

% nSequence: [working image size, max iteration at this size]
%ri.nFinal=72; %%%%%%
ri.nSequence=[32 7
             48 12
             72 20
      ri.nFinal ri.nIters];
ri.nSequence=[48 ri.nIters];
% ri.nSequence=[48 5
%               72 10
%        ri.nFinal ri.nIters]; %%%%%%%

ri.mergeIter=inf;  % First iteration at which twin volumes are merged
ri.highResIter=1; % k is higher before this iteration
ri.altImgIter=ri.nIters-1;
% ri.altImgIter=inf;
ri.resetActiveRVIters=[14 21];  % Force a full search at this iteration
ri.skewIter=0;

% Assign the image size according to the start iteration number.
iseq=find(ri.nSequence(:,2)>=ri.startIter,1);
        if numel(iseq)<1
            iseq=1;
        end;
ri.nCurrent=ri.nSequence(iseq,1);


% Partitioning the reconstruction among nodes and workers
if ri.nVols>1
    ri.nTwins=1;   % Do both odd and even
else
    ri.nTwins=normalNTwins;
end;
    % ri.twinMode='RsoIso';
ri.twinMode='OddEven';
ri.nGroups=1;  % default, to be updated by number of jobs
ri.nSlices=1;  % number of slices per job

%------------special cases for 2D classification-----------
if mode2D
    ri.altImageIter=inf;
    ri.nVols=1;
    ri.nTwins=1;
end;


% Iterations and filenames
% if isShort
%     ri.nIters=7;
% end;

%  It's at the last iteration that the altImgs (membrane-restored) are used

% Waiting times: 1. for moi (initial sync) 2. for roi, fvs (output of EM)
%                3. for fsc (sync for twin reconstructions.
ri.baseTimeout=[30 30 30]*60*ri.nVols; % waiting times in s at 64 pixel size

if isSmall
    ri.baseTimeout=ri.baseTimeout/5;
end;
    if isShort
        ri.baseTimeout=[30 30 30]*10*ri.nVols;
    end;

% Stack file to use
ri.stackPath='Stack/';
% ri.stackPath='StackMerge/';
ri.siName='sq03_1fp128tsi.mat'; % 151117
ri.activeFlagSet=10;  % strict manual picking, 13k particles
%ri.activeFlagSet=9;  % vol1 from subsequent 72k5v4

% ri.siName='Merge1415tsi.mat'; % 140625 + 151117 merged stack
% ri.siName='sq10_350w10-00p128vtsi.mat'  % First exposure only
% ri.siName='sq04_1288p128v11tsi.mat';
% ri.siName='114119asi.mat';
%ri.siName='Simn1000siz128sigma99clk00tsi.mat';
%ri.activeFlagSet=1;  % vol1 from subsequent 72k5v4

% Output files
pathBasename=[num2str(ri.nFinal) pathSuffix];
if mode2D
    ri.outPath=['Reconstructions/Class' pathBasename '/']; % relative to the root path
else
    ri.outPath=['Reconstructions/Recon' pathBasename '/']; % relative to the root path
end;
ri.tempPath=[ri.outPath 'temp/'];
ri.outBasename='';  % Prefix for output files

% Reconstruction Wiener constant
ri.kFactor=4;

if isSmall
    ri.kFactor=1;
end;

% ------Load the raw stack info and set the ri.activeFlags-------
si=LoadStruct([ri.stackPath ri.siName]);
nSets=size(si.activeFlags,2);
nParts=size(si.miIndex,1);

if size(si.activeFlags,1) ~= nParts
    error('Inconsistent size of si.activeFlags');
end;

% form the union of the activeFlags
ri.activeFlags=false(nParts,1);
for i=1:numel(ri.activeFlagSet)
    if ri.activeFlagSet(i)>nSets
        warning(['ri.activeFlagSet too large: ' num2str(ri.activeFlagSet(i)) ...
            ' changed to ' num2str(nSets)]);
        ri.activeFlagSet(i)=nSets;
    end;
    ri.activeFlags=ri.activeFlags | si.activeFlags(:,ri.activeFlagSet(i));
end;

if isShort
    ri.activeFlags(2001:end)=false;
end;

% Set up the twin flags
ri.twinFlags=false(nParts,ri.nTwins);
switch lower(ri.twinMode)
    case 'oddeven'
        activeInds=find(ri.activeFlags);
        for i=1:ri.nTwins
            twinInds=activeInds(i:ri.nTwins:end);
            ri.nTwin(1,i)=numel(twinInds);
            ri.twinFlags(twinInds,i)=true;
        end;
    case 'rsoiso'
        %         Find all the rso particles
        rsoFlags=false(nParts,1);
        for j=1:numel(si.miIndex)
            imi=si.miIndex(j);
            ip=si.miParticle(j);
            rsoFlags(j)=si.mi{imi}.particle.picks(ip,7);
        end;
        ri.twinFlags(:,1)=rsoFlags & ri.activeFlags;
        ri.twinFlags(:,2)=(~rsoFlags) & ri.activeFlags;
    otherwise
        error(['Unrecognized twinMode: ' ri.twinMode]);
end;

ri.nTwin=sum(ri.twinFlags,1);  % Number of images in each twin set.

% % defaultMbnOffsetA=-70;  % position of membrane relative to RSO particle center
    % 3D starting models
    % Initial volumes are put into moi.
    % First, figure out our final pixel size from the si structure.
    ri.pixAU=si.pixA;  % pixel size of original stack
    pixA=ri.pixAU*ri.nCropU/ri.nCurrent;
    % Get the moi.refVols starting volumes
    [origVols,ri.mbnOffsetA]=arGetRefVolumes(pixA,ri.nCurrent,ri.startMapName,ri.nVols);
    moi.refVols=zeros(size(origVols),'single');
    fc=pixA/ri.startRes;  % lowpass for starting volume
    ri.volSD=zeros(1,ri.nVols);
    for iVol=1:ri.nVols
        v=SharpFilt(origVols(:,:,:,iVol),fc,fc/ri.nCurrent);
        ri.volSD(iVol)=sqrt(v(:)'*v(:)/numel(v));  % normalize the volume
        %     moi.refVols(:,:,:,iVol)=v/ri.volSD(iVol);
        moi.refVols(:,:,:,iVol)=ri.refScale*v;
    end;

if mode2D
    ri=reSetRefAngles(ri.angleStepsU,ri.angleLimits,ri.isos,false,ri);
    ri.angles=reGetAngleList(ri,false); % Count the references    
    n=ri.nCurrent;
    nRefs=size(ri.angles,1);
    
    moi.refs=zeros(n,n,nRefs,'single');
    fcr=fc;
    if ri.flags.useRandomRefs<2
        fcr=fc*2;  % higher-frequency noise if adding to blob refs
    end;
    if ri.flags.useRandomRefs % nonzero flag
        refMask=repmat(fuzzymask(n,2,0.3*n,0.1*n),1,1,nRefs);
        moi.refs=SharpFilt(randn(n,n,nRefs).*refMask,fcr);
        moi.refs=moi.refs/std(moi.refs(:));
    end;
    
    if ri.flags.useRandomRefs<2  % flag is zero or 1: use projections of model
        % Initial volumes and are put into moi; their projections are the initial
        % values of the classes.
        % First, figure out our final pixel size from the si structure.
        % Make the 2D refs
        moi.refs=moi.refs+reMakeTemplates(moi.refVols,ri.angles);  % refs(x,y,iRef,iVol)
    end;
end;

