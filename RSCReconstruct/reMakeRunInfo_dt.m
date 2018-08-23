function [ri,si,moi]=reMakeRunInfo
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

pathSuffix='dt';

% ---- Flags ----
fl=struct;
fl.makeFakeData=false;
fl.doShowPriors=false;
fl.makeRawProbs=false;
fl.useAltImgs=true;  % use the unsubtracted stack for alternative reconstruction.
fl.useAltImgs=~isShort; %%% dont use them for short problem.
% fl.removeRings=true;
fl.removeRings=false;


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
moi.imgAmp=1;
moi.b0=0;
moi.pRefs=1;
moi.refVols=0;  % placeholder
moi.pVols=1;    % placeholder
moi.intensiveFields=fieldnames(moi);

% ----------Create the run info----------
ri=struct;
ri.flags=fl;
ri.usePrior=1;
ri.startMapName='KvMap.mat';
ri.startRes=60;

ri.thresholds=[1e-5 0 0];  % Thresh for trans, alpha, refs
ri.refScale=.1;     % reference amplitude factor
ri.symmetry=4;

% Angle and shift search parameters, and start assigning the run info ri
ri.isos=[0 1]; % to handle rso and iso particles.
ri.maxShiftU=6; % Maximum translation searched

ri.maxTxFraction=.7;
ri.angleStepsU=[3 3 5];  % alpha, beta, gamma steps in degrees

% if isSmall || isShort
%    ri.maxShift=4;
%    ri.angleSteps=[6 6 10];
% end;

ri.angleLimits=[20 0 360/ri.symmetry];  % +-alpha range;  beta start; max gamma.

% Number of volumes
ri.nVols=1;
moi.pVols=ones(ri.nVols,1);

ri.startIter=1;
ri.nIters=31;

% Final image sizes
ri.nCropU=112;  % We first crop the input images to this size
ri.nFinal=ri.nCropU;

% nSequence: [working image size, max iteration at this size]
ri.nSequence=[32 6
              48 11
              72 20
       ri.nFinal ri.nIters];

ri.mergeIter=inf;  % First iteration at which twin volumes are merged
ri.highResIter=10; % k is higher before this iteration
ri.altImgIter=ri.nIters-1;  %%%% start this at 30

% Assign the image size according to the start iteration number.
iseq=find(ri.nSequence(:,2)>=ri.startIter,1);
        if numel(iseq)<1
            iseq=1;
        end;
ri.nCurrent=ri.nSequence(iseq,1);


% Partitioning the reconstruction among nodes and workers
ri.nTwins=2;   % Do both odd and even
% ri.twinMode='RsoIso';
ri.twinMode='OddEven';
ri.nGroups=1;  % default, to be updated by number of jobs
ri.nSlices=1;  % number of slices per job

% Iterations and filenames
if isShort
    ri.nIters=7;
end;

%  It's at the last iteration that the altImgs (membrane-restored) are used

% Waiting times: 1. for moi (initial sync) 2. for roi, fvs (output of EM)
%                3. for fsc (sync for twin reconstructions.
ri.baseTimeout=[30 30 30]*60; % waiting times in s at 64 pixel size

if isSmall
    ri.baseTimeout=ri.baseTimeout/5;
end;
    if isShort
        ri.baseTimeout=[100 100 100];
    end;

ri.activeFlagSet=4;  % In this stack, this second pruning

% Stack file to use
ri.stackPath='Stack/';
ri.siName='sq03_1fp128tsi.mat';
% ri.siName='sq10_350circfp128tsi.mat';
% ri.siName='sq10_350nosubtsi.mat'; ri.activeFlagSet=4;
% ri.siName='sq10_350p128tsi.mat'; ri.activeFlagSet=4; original 140625 stack
% ri.siName='sq10_350w10-00p128vtsi.mat'  % First exposure only
% ri.siName='sq04_1288p128v11tsi.mat';
% ri.siName='114119asi.mat';

% Output files
pathBasename=[num2str(ri.nFinal) pathSuffix];
ri.outPath=['Reconstructions/Recon' pathBasename '/']; % relative to the root path
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
if ri.activeFlagSet>nSets
    warning(['ri.activeFlagSet too large: ' num2str(ri.activeFlagSet) ...
        ' changed to ' num2str(nSets)]);
    ri.activeFlagSet=nSets;
end;
ri.activeFlags=si.activeFlags(:,ri.activeFlagSet);

if isShort
    ri.activeFlags(201:end)=false;
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


