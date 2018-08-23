function [templateAngles ri]=rsMakeTemplateAngles3(ri)
% Given a RunInfo structure, make a list of beta and gamma angles.
% return the nangs x 3 array, with gamma varying most quickly.
% The ri.gammaMode selects various ways to handle the density of points
% near the poles of the sphere.  gammaMode='uniform' ignores the problem
% and leaves the number of gamma values independent of beta.
% gammaMode=variable has the number of gamma values follow sin(beta).
% gammaMode=quantized is the same, except that the number of gamma values
% changes in such a way that no new gamma values appear in the process.
% This is useful if rotations through gamma are expensive.
% 
%   For example, if the following fields of ri are filled,
% ri.nSteps=[40 20 20]; % no. of steps in alpha, beta and gamma
% ri.symmetry=   2;  % We're assuming C symmetry only
% ri.gammaMode='quantized';
% 
%   the returned copy of ri will also have these fields filled:
% ri.angleSteps=[0 9.4737 0]  % in degres
% ri.nGammas (20 x 1)  % no. of gammas at each beta step
% ri.dGammas [20 x 1]  % corresponding gamma step
% ri.iBetalookup=[370 x 1]; % beta index at each angle index
% ri.angleIndexLookup= [20 x 1] % first angle index for the given beta
% index.
% 
if ~isfield(ri,'gammaMode')
    ri.gammaMode='uniform';
end;

nBetas=ri.nSteps(2);
ri.angleSteps(1)=360/ri.nSteps(1);
ri.angleSteps(2)=180/(nBetas-1);  % correct the step size.
iBetas=(1:nBetas)';
betas=(iBetas-1)*ri.angleSteps(2);

nGamma0=ri.nSteps(3);
ri.angleSteps(3)=360/(ri.symmetry*nGamma0);  % correct the step size
switch lower(ri.gammaMode)
    case 'uniform'
        ri.nGammas=repmat(nGamma0,nBetas,1);
    case 'variable'
        ri.nGammas=max(1,ceil(sind(betas)*nGamma0));
    case 'quantized'
        for i=iBetas'
            stride=floor(min(nGamma0,1/sind(betas(i))));
            %         force the stride to be a factor of nGamma0
            while (stride>1) && (mod(nGamma0,stride)~=0)
                stride=stride-1;
            end;
            ri.nGammas(i)=nGamma0/stride;
        end;
    otherwise
        error('unrecognized gammaMode');
end;
ri.dGammas=360./(ri.nGammas*ri.symmetry);

nangs=sum(ri.nGammas);  % total nuymber of entries
ri.iBetaLookup=single(zeros(nangs,1));
ri.angleIndexLookup=single(zeros(nBetas,1));
%%
templateAngles=zeros(nangs,3);
i=1;
for iBeta=1:nBetas
    ri.angleIndexLookup(iBeta)=i;
    i0=i;
    for iGamma=1:ri.nGammas(iBeta)
        templateAngles(i,:)=[0  betas(iBeta) (iGamma-1)*ri.dGammas(iBeta)];
        i=i+1;
    end;
    ri.iBetaLookup(i0:i-1)=iBeta;
end;
