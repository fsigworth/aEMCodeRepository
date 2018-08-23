function [refAngles, isoList]=reGetAngleList(ri,makeAlphas)
% function [refAngles, isoList]=reGetAngleList(ri,makeAlphas)
% Given the ri structure, get the list of reference angles.  The isoList
% has the same number of rows as refAngles, but is either 0 or 1 according
% to which orientation the particle will have.
if nargin<2
    makeAlphas=0;
end;
if makeAlphas  % set up the working alphas list.
    alphasW=ri.alphas(:);  % This must have twice as many elements as ri.alphas
    nAlphas=size(ri.alphas,1);
else
    alphasW=0;
    nAlphas=1;
end;
    nBetaGammas=prod(ri.angleN(2:3));
    [betas, gammas]=reGetBetaGammaList(ri);
    nAlphasW=numel(alphasW);
    alphas2=repmat(alphasW(:),1,nBetaGammas);
    betas2=repmat(betas,nAlphasW,1);
    gammas2=repmat(gammas,nAlphasW,1);
    refAngles=[alphas2(:) betas2(:) gammas2(:)];
%     isoList=repmat(isos(:)',numel(ri.alphas),nBetaGammas,numel(ri.alphas));
if nargout>1  % create the list of angles, including isos
    isoList=repmat(ri.isos(:)',nAlphas,nBetaGammas);
%     if ri.isos has 2 elements, isoList is nAlphas x 2 x nBetaGammas in
%     size.
    isoList=isoList(:);
end;
