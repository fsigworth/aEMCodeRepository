% RSCDefs2
% Definitions of RSC datastructures and functions
% v2.0  6 Aug 2012

%% DEFINE THE DATASET
% We assume that we have loaded into memory an array of micrographInfo structures according to
function [miArray basenames]=rsLoadMicrographInfos(miPath,basenames)
        % load every mi file in the given path.  If basenames (cell array of strings)
        % is given, load the files according to the sequence of basenames.
        % Or if desired the list of basenames is returned.

%% Description of the run: iterations with the same set of basic parameters
runInfo
    ri.pixA          =   3   % angstroms per pixel/voxel in reconstruction, may differ from micrographs
    ri.n             =  64  % box size
    ri.membraneOffset= -13  % z-coordinate of membrane center in the map.  Negative for AMPAR, positive for BK.
    ri.maxTranslation=   7 % size of the cross-correlation search is 2*maxTranslation+1
    ri.angleSteps     =  [5 5 5]  % angular steps in degrees.
    ri.symmetry      =   2  % We're assuming C symmetry only
    ri.alphaMode = 'uniform';
    ri.gammaMode = 'variable';   % options: quantized, variable
    % The following fields are created by rsMakeTemplateAngles3:
    ri.ibetaLookup(angIndex)    % Given index into the angles, get the ibeta
    ri.betaAngleIndex(ibeta)    % Given ibeta, get the first index into the angles
    ri.nAlphas(ibeta);
    ri.nGammas(ibeta)           % Number of gamma values for this beta

%     
% % templateAngles are stored in an nang x 3 array, where all alphas are zero,
% %     betas change most slowly, and gammas change quickly.  To get the beta
% %     and gamma values corresponding to angleIndex do the following:
%     ibeta=ai.ibetaLookup(angleIndex);
%     beta=ibeta*ri.angleSteps(2);
%     gammaStep=360/(ri.symmetry*ai.nGammas(ibeta));
%     indexOffset=angleIndex-ai.betaAngleIndex(ibeta);
%     gamma=indexOffset*gammaStep;
%     
%     switch ai.gammaMode
%         case 'uniform'
%             igamma=diffIndex;
%         case 'quantized'
%             igamma=diffIndex * ai.gammaStep;
%         case 'variable'
%             igamma=diffIndex
% function [ibeta igamma]=rsInd2Sub(angleIndex)
% function [angleIndex ngamma]=rsSub2Ind(ibeta)

%% TEMPLATES
% Angles are defined as [alpha beta gamma].  Looking down on the north pole,
% where alpha is longitude (prime meridian is y axis), increasing cw; beta is latitude 
% (0 = north pole) and gamma is local self-rotation (cw viewed from top of particle,
% when particle is at [0 0 0]).
% Orientation is defined as [alpha beta gamma iso lower] where iso=1 for
% inside-out particle, lower=1 for lower hemisphere.
% note that [alpha beta gamma 1 0]=[alpha+180, 180-beta gamma 0 0] with a shift of magnitude 2b*sin(beta).
%           [alpha beta gamma 0 1]=[alpha 180-beta gamma 0 0]
%           [alpha beta gamma 1 1]=[alpha+180 beta gamma 0 0] with a shift.
% Functions operating on angles can take m x 3 matrices as vectors of angles.
% Templates are stored simply as a stack, the array n x n x nTemplates, with no CTF filtering.
function eulerAngles= rsDegToEuler(rscAngles)
% Converts [alpha beta gamma] in degrees to Euler angles in radians
% 
% function templateAngles=rsMakeTemplateAngles(runInfo)  % return rscAngles, with gamma always 0.
function [templateAngles ri]=rsMakeTemplateAngles3(ri)  % return rscAngles (nangs x 3) and fill in the 
% runInfo fields for looking up angles.  templateAngles has all values of
% beta and gamma, but alpha is returned as zero.
function templates=rsMakeTemplates(templateAngles,map)  % templates is array n x n x nAngles
% function imageLoc=rsParticlePosition(templateAngles)   % for each viewing direction, compute
% templateAngles is a matrix nTemplates x 3, of angles in degrees [alpha beta 0]
% function tIndices=getNearestTemplate(runInfo, templateAngles, angles)


%% PARTICLE STACK
% From the micrograph info structures and micrographs that we read in,
% we make a stack of particle images and
% a stackInfo structure that contains all the relevant information for reconstruction, plus
% pointers back to the original image.

stackInfo
    si.pixA         % Copied from runInfo
    si.miIndex(i)   % Index of the micrograph corresponding to particle i (index of the micrographInfos)
    si.miParticle(i)% Index of the particle within the micrograph.  This allows lookup of x,y,vesR etc.
    si.alpha0(i)    % alpha rotation applied to the image
    si.yClick(i)    % Copied from the mi struct.  xClick is zero after rotation.
    si.rVesicle(i)  % Copied from the mi struct

function [stackInfo stack]=rsMakeParticleStack(runInfo,miArray,startIndex,endIndex,type,quality)
    % This function boxes (and resamples to the correct pixel size if necessary)
    % the particles from the original images.
    % The optional startIndex and endIndex allows a subset of the stack to be generated.
    % Optionally, only particles that match the type and quality values stored in the mi structures
    % are stored.
    % This function applies the alpha rotation computed from the particle position.

function [ctfs]=rsMakeCTFStack(miArray,stackInfo)
    % This function produces a stack of ctfs, one ctf per relevant micrograph.
    % The ctf for particle i will be given by ctfs(:,:,stackInfo.miIndex(i)).

%% UPDATES FOR EVERY ITERATION

% 1. Best-match variables per image.  This struct can also define a simulated data set.
%   This can be updated in each iteration, or used for best-match reconstruction.
particleInfo
    pa.bestRefIndex(i)
    pa.bestGammaIndex(i)
    pa.translation(i)   %(2x1 vector, in pixels)
    pa.alpha(i)         %(additional alpha rotation)
    pa.amp(i)           %(signal amplitude) 
    pa.active(i)        %(0 or 1)

% 2.  Global model parameters, which along with the map are updated with
% each EM iteration.  This everything in X except for the 3D volumes

modelInfo
    mo.noiseSD    % noise standard deviation
    mo.clickSD    % in pixels
    mo.rockSD     % composite rocking and bobbing
    mo.bobSD
    mo.bobMean
    mo.pIO        % probability inside-out
    mo.pVols(k)      % probability of the various 3D volumes
    mo.likelihood

%% DOMAIN AND PROBABILITIES

% Based on the prior, click error etc. we want to define an integration domain for each image.  Set a threshold
% such that values of the prior that fall below that value are not considered.  We then select a subset
% of templates to use.  We also select a subset of psi angles.

domain
    .templateList   % nTL x 1 array containing tIndex values of the templates we'll use
    .alphaList      % nPL x 1 array containing integral alpha values that we'll use.
    .nCC            % size of the cross correlations we preserve
    .indices(j,:)=[jTL,jPsiL,jIO]  % Description of the domain.  This contains 
                    % the indices into the lists for each term that we'll use.
        % psiList is a list of nPI iPsi values, where the psi angles through which the image is to be
        % rotated are psi=iPsi*thetaStep.  This allows pre-computation of rotated images.

function [domain pPriorAndClick]=GetDomainFromClick(mi,stackInfo,runInfo,particleInfo,i,domainThreshold)
    % Find the domain for computing data probabilities for image i, based on the
    % prior and click probabilities.
    % We keep only terms where their product is above the domainThreshold.  pPriorAndClick(:,:,j) contains
    % the corresponding nCC x nCC translation probabilities for each template, psi and iO
    % combination.

function [pLatent aShift aVar]=GetLatentProbability(image,ctf,templates,pPriorAndClick,domain)
    % Given the particle image and its ctf, compute the latent probability function
    % and also the expectation values aShift for the estimation of the b0 displacement, and aVar for
    % re-estimation of the noise variance.
    % The latent probability is an array of size (domain.nCC^2 x size(domain.indices,1)), same as
    % pPriorAndClick.

function [classIncrements classNorms]=GetClassIncrements(image,ctf,model,pLatent)
    % Get a singe term for the sums D and zeta (eqns. 2.18 and 2.19) from the image given.
    %  Each is n x n x nTemplates in size.

