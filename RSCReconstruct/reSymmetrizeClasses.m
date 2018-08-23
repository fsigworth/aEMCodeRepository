function out=reSymmetrizeClasses(in,ri)
% For ri.symmetry is even, average over the symmetry such that projections
% at (alpha 180-beta gamma) are forced to be mirror images of those at
% (alpha beta gamma).
% in: class means.  They are assumed to n x n x nRefs x nVols x nTwins
% ri: the run info structure
% The ref angles all have alpha=0 and gamma changes faster than beta.


if mod(ri.symmetry,2)==0 % Support only even symmetries
    % Enforce the fact that [alpha beta gamma] is mirror image of
    %[-alpha 180-beta gamma].  Here we assume alpha is always zero.
    out=in;
    nBeta=ri.angleN(2);
    nGamma=ri.angleN(3);
    for iBeta=1:nBeta
        imgs=in(:,:,(iBeta-1)*nGamma+1:iBeta*nGamma,:,:,:)...
            +MirrorX(in(:,:,(nBeta-iBeta)*nGamma+1:(nBeta-iBeta+1)*nGamma,:,:,:));
        out(:,:,(iBeta-1)*nGamma+1:iBeta*nGamma,:,:,:)=imgs/2;
        out(:,:,(nBeta-iBeta)*nGamma+1:(nBeta-iBeta+1)*nGamma,:,:,:)=MirrorX(imgs/2);
    end;
else
    error(['This symmetry not supported: ' num2str(ri.symmetry)]);
end;
