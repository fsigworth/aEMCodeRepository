function [ri,si,moi]=re2DMakeRunInfo
% Sets up the structures ri, si and moi.
% ri is the run info for all workers; si is the entire raw stack; moi is
% the initial model. We assume that we in a directory containing
% Reconstructions/ and Stack/ .  Within Reconstructions/ the main program
% will create directories such as 96a/ where 'a' is given by pathSuffix.

% isSmall=false;  % run a small problem, e.g. on laptop
% isSmall=true;
isSmall=false;
% isShort=true;
isShort=false;

pathSuffix='s10c1r1';
mode2D=false;

% ---- Flags ----
fl=struct;
fl.mode2D=mode2D;
fl.makeFakeData=false;
fl.doShowPriors=false;
fl.makeRawProbs=false;
fl.useAltImgs=true;  % use the unsubtracted stack for alternative reconstruction.
fl.useAltImgs=~isShort; %%% dont use them for short problem.
fl.useRandomRefs=2;  % 0: model projections; 1: model+random; 2: random only

fl.writeAllGroupFiles=false;
fl.saveRoi=true;  % always save the composite roi structure for each twin
fl.saveMoi=true;   % save the composite moi structures
fl.saveFigures=true;
fl.useParFor=false;

% ---------Create the initial model---------
moi=struct;  % initialize the intensive model variables
moi.sigmaN=10;  % start with a high value
moi.sigmaC=2;
moi.sigmaG=2;
moi.imgAmp=1; % starting value
moi.b0=0;
moi.pRefs=1;
moi.refVols=0; % placeholder
moi.refs=0;  % placeholder
moi.intensiveFields=fieldnames(moi);

% ----------Create the run info----------
ri=struct;
ri.flags=fl;

ri.startMapName='KvMap.mat';
ri.startRes=60;

ri.thresholds=[1e-5 0 0];  % Thresh for trans, alpha, refs
ri.refScale=.1;     % reference amplitude factor
ri.symmetry=4;
ri.usePrior=1;

% Angle and shift search parameters, and start assigning the run info ri
ri.isos=[1 1]; % to handle rso and iso particles.
ri.maxShiftU=6; % Maximum translation searched
ri.maxTxFraction=.7;
ri.angleStepsU=[5 5 8];  % alpha, beta, gamma steps in degrees
if isShort
    ri.angleStepsU=[10 20 20];
end;
% if isSmall || isShort
%    ri.maxShift=4;
%    ri.angleSteps=[6 6 10];
% end;

ri.angleLimits=[10 0 360/ri.symmetry];  % +-alpha range;  beta start; max gamma.

% Number of volumes
ri.nVols=1;
moi.pVols=ones(ri.nVols,1);

ri.startIter=1;
ri.nIters=30;

% Working image sizes
ri.nCropU=108;  % We first crop the input images to this size
% ri.nFinal=96;
ri.nFinal=64;

% nSequence: [working image size, max iteration at this size]
% ri.nSequence=[32 6
%               48 11
%               64 24
%        ri.nFinal ri.nIters];

ri.nSequence=[ ri.nFinal 24
    ri.nFinal ri.nIters];


ri.mergeIter=20;  % First iteration at which twin volumes are merged
ri.highResIter=10; % k is higher before this iteration
ri.altImgIter=ri.nIters;

ri.nCurrent=ri.nSequence(1,1);

% if isSmall
%     ri.nU=64;
%     ri.nCrop=48;
% end;

% Partitioning the reconstruction among nodes and workers
 ri.nTwins=1;   % No twins for 2D classification
% ri.twinMode='RsoIso';
ri.twinMode='OddEven';
ri.nGroups=2;  % default, to be updated by number of jobs
ri.nSlices=1;  % number of slices per job

%  It's at the last iteration that the altImgs (membrane-restored) are used

% Waiting times: 1. for moi (initial sync) 2. for roi, fvs (output of EM)
%                3. for fsc (sync for twin reconstructions.
ri.timeout=[10 30 30]*60; % waiting times in s

% if isSmall
%     ri.timeout=ri.timeout/5;
% end;
%     if isShort
%         ri.timeout=[1 1 1];
%     end;

% Stack file to use
ri.stackPath='Stack/';
ri.siName='sq03_1fp128tsi.mat';

% % ri.siName='sq10_350mlfp128tsi.mat';
% % ri.siName='sq10_350p128tsi.mat'; % original 140625 stack
% ri.siName='sq10_350circfp128tsi.mat'; % original 140625 stack
% % ri.siName='sq10_350w10-00p128vtsi.mat'  % First exposure only
% % ri.siName='sq04_1288p128v11tsi.mat';
% % ri.siName='114119asi.mat';

% ri.activeFlagSet=2;  % In this stack, this is both iso and rso particles.
ri.activeFlagSet=10;  % 

% Output files
pathBasename=[num2str(ri.nFinal) pathSuffix];
ri.outPath=['Reconstructions/Class' pathBasename '/']; % relative to the root path
ri.tempPath=[ri.outPath 'temp/'];
ri.outBasename='';  % Prefix for output files

% Reconstruction Wiener constant
ri.kFactor=4;

% if isSmall
%     ri.kFactor=1;
% end;

% ------Load the raw stack info and set the ri.activeFlags-------
si=LoadStruct([ri.stackPath ri.siName]);
nSets=size(si.activeFlags,2);
nParts=size(si.miIndex,1);

if size(si.activeFlags,1) ~= nParts
    error('Inconsistent size of si.activeFlags');
end;
if ri.activeFlagSet>nSets
    warning(['ri.activeFlagSet too large: ' num2str(ri.activeFlagSet) ...
        ' changed to ' num2str(nSets)]);
    ri.activeFlagSet=nSets;
end;
ri.activeFlags=si.activeFlags(:,ri.activeFlagSet);

if isShort
    ri.activeFlags(3001:end)=false;
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


ri=reSetRefAngles(ri.angleStepsU,ri.angleLimits,ri.isos,false,ri);
ri.angles=reGetAngleList(ri,false); % Count the references
ri.pixAU=si.pixA;  % pixel size of original stack
pixA=ri.pixAU*ri.nCropU/ri.nCurrent;
fc=pixA/ri.startRes;  % lowpass for starting volume

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
    % Get the moi.refVols starting volumes
    [origVols,ri.mbnOffsetA]=arGetRefVolumes(pixA,ri.nCurrent,ri.startMapName,ri.nVols);
    moi.refVols=zeros(size(origVols),'single');
    
    ri.volSD=zeros(1,ri.nVols);
    for iVol=1:ri.nVols
        v=SharpFilt(origVols(:,:,:,iVol),fc,fc/ri.nCurrent);
        moi.refVols(:,:,:,iVol)=ri.refScale*v;
    end;
    
    % Make the 2D refs
    moi.refs=moi.refs+reMakeTemplates(moi.refVols,ri.angles);  % refs(x,y,iRef,iVol)
end;
end



