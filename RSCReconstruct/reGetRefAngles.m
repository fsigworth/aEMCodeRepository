function [betas, gammas]=reGetRefAngles(ri,refInds)
% function [betas, gammas]=reGetRefAngles(ri,refInds)
% *** Deprecated function, use reGetBetaGammaList ***
% From the ri structure, give a list of the betas and gammas.  If the
% optional refInds is given, give only the angles corresponding to those
% indices.
% The references are enumerated with gamma varying most quickly.
% The ri structure contains the following 3-element vectors, referring to
% alpha, beta, gamma (in degrees)
% ri.angleMins
% ri.angleSteps
% ri.angleN
% 
maxInds=prod(ri.angleN(2:3));
if nargin<2
    refInds=1:maxInds;
end;
[betas, gammas]=reGetBetaGammaList(ri,refInds);
