function [ri, refAngles, isoList]=reSetRefAngles(angleSteps,angleLimits,isos,makeAlphas,ri)
% function [ri,refAngles,isoList]=reSetRefAngles(angleSteps,angleLimits,isos,makeAlphas,ri)
% angleLimits is defined as [minAlpha, minBeta, maxGamma].
% e.g. angleSteps=[10 15 20], angleLimits=[-20 30 180]
% will result in 5 alphas (always odd) from -20 to 20,
% 8 betas (always even)from 30 to 150 (with 17.14 degree steps),
% and 9 gammas from 0 to 160.
%
% To get only 1 beta equal to angleLimits(2),
%   set angleSteps(2)>2*(90-angleLimits(2)).
% To get only one gamma=0, set angleLimits(3)=0, angleSteps(3) nonzero.
%
% isos=0 -> only rso particles will be considered; isos=1 -> only iso
% particles; isos=[0 1] (default) -> both, and twice as many alphas will be
% generated.
% If no input ri is given, ri is created.  It has these 3-element vectors:
% ri.angleMins
% ri.angleSteps
% ri.angleN
% plus the 1 or 2 element vector ri.isos copied from the input argument.
%
% If makeAlphas==0 (default) the refAngles=[0 betas gammas] will be
% created, where nRefs=ri.angleN(3)*ri.angleN(2), and gammas change most
% rapidly.
% If makeAlphas is set, then alphas vary too,
% more rapidly than isos, then gammas, then betas change most slowly.
% All angles are created according to the following (i is angle index, j is
% which angle):
%   angle=ri.angleMin(j)+(i-1)*ri.angleStep(j), with iso*180 added to alpha.
%  The returned refAngles are ready to be used in rsMakeTemplates.    
if nargin<3
    isos=[0 1];  % Handle rso and iso
end;
if nargin<4
    makeAlphas=0;
end;

if nargin<5
    ri=struct;
    ri.isos=isos;
end;

% Construct the list of alphas
% There will be an odd number of alphas, including zero.
ni=numel(isos);
ri.isos=isos;
if abs(angleLimits(1))>0 && abs(angleSteps(1)>0)
    nAlpha2=round(abs(angleLimits(1)/angleSteps(1)));
    ri.angleStep(1)=abs(angleLimits(1))/nAlpha2;
    ri.angleMin(1)=-abs(angleLimits(1));
    nAlphas=2*nAlpha2+1;
    ri.angleN(1)=nAlphas;
    alphas=ri.angleMin(1)+(0:nAlphas-1)*ri.angleStep(1);
    ri.alphas=zeros(nAlphas,ni);
    for i=1:ni
        ri.alphas(:,i)=alphas+isos(i)*180;
    end;
else
    ri.angleMin(1)=0;
    ri.angleStep(1)=0;
    ri.angleN(1)=1;
    ri.alphas=zeros(1,ni);
end;

% There will be an even no. of betas, bracketing 90 degrees.
% To get only 1 beta equal to angleLimits(2),
%   set angleSteps(2)>2*(90-angleLimits(2)).
ri.angleMin(2)=angleLimits(2);
nBeta2=round((90-angleLimits(2))/angleSteps(2));
ri.angleStep(2)=(90-angleLimits(2))/(nBeta2-0.5);
ri.angleN(2)=max(1,2*nBeta2);  % always give one beta

% The angleLimits(3) gives the max. gamma, typically 360/symmetry.
% To give only one gamma=0, set angleLimits(3)=0, angleSteps(3) nonzero.
ri.angleN(3)=max(1,round(angleLimits(3)/angleSteps(3)));  % at least 1.
ri.angleStep(3)=angleLimits(3)/ri.angleN(3);
ri.angleMin(3)=0;

if nargout>1
    [refAngles, isoList]=reGetAngleList(ri,makeAlphas);
end;
