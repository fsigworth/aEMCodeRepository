function [out,reconAngles]=reSymmetrizeClassesX(in,ri)
% in: class means.  They are assumed to n x n x nRefs x nVols x nTwins
% ri: the run info structure
% Converts the class means in into projections for reconstruction, modifying
% for cyclic symmetries, and generating an expanded list of angles for
% reconstruction, in which gammas are incremented.  
% nRefs is assumed equal to ri.angleN(2)*ri.angleN(3).
% For even ri.symmetry, average over the symmetry such that projections at
% (alpha 180-beta gamma) are forced to be mirror images of those at (alpha
% beta gamma).  For higher symmetries the output is increased in size to
% n x n x nRefs*symmetry/2 x nVols x nTwins


riSym=ri;
[betas,gammas]=reGetRefAngles(ri);
alphas=0*betas;
reconAngles=[alphas betas gammas];
if ri.symmetry>1
    % Enforce the fact that [alpha beta gamma] is mirror image of
    %                                   [-alpha 180-beta gamma].  Here we
    %                                   assume alpha is always zero.
    nBeta=ri.angleN(2);
    nGammas=ri.angleN(3);
    for iBeta=1:nBeta
        imgs=in(:,:,(iBeta-1)*nGammas+1:iBeta*nGammas,:,:,:)...
            +MirrorX(in(:,:,(nBeta-iBeta)*nGammas+1:(nBeta-iBeta+1)*nGammas,:,:,:));
        out(:,:,(iBeta-1)*nGammas+1:iBeta*nGammas,:,:,:)=imgs/2;
        out(:,:,(nBeta-iBeta)*nGammas+1:(nBeta-iBeta+1)*nGammas,:,:,:)=MirrorX(imgs/2);
    end;
else
    out=in;
    return
end;

if mod(ri.symmetry,2)==0  % Even symmetry; need to extend the gamma values
    sym2=ri.symmetry/2;
    gammaLimit=gammaStep*nGammas;
    if gammaLimit>360/symmetry  % already sampling too many
        warning(['Symmetry ' num2str(ri.symmetry) ' but gamma sampled up to ' num2str(gammaLimit)]);
        return
    end;
    riSym.angleN(3)=sym2*ri.angleN(3);
    reps=ones(1,ndims(in));
    reps(3)=sym2;  % vector of replicated dimensions
    out=repmat(in,reps);
    for i=1:sym2-1
        newGammas=gammas+i*180/sym2;
        reconAngles=[reconAngles; alphas betas newGammas];
    end;
else
    error(['This symmetry not supported: ' num2str(ri.symmetry)]);
end;
