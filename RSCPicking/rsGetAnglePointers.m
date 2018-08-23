function ptrs=rsGetAnglePointers(angles,nAlpha,nBeta,listIndices)
% Given a t x 2 array angles, find a t x 1 array of pointers.
% Each pointer gives the closest index of the angleList created by
% rsListHemisphereAngles, using its listIndices lookup table.
%
corr=(mod(nBeta,2)+1)/2; % correction is .5 even, 1 odd
betaStep=90/(nBeta-corr); % odd nBeta gives a final value at 90 degrees;
% even nBeta gives a final value of 90*(nBeta-1)/(nBeta-.5))
alphas=mod(angles(:,1),360);
betas=mod(angles(:,2),180);
% while any(alphas<0)
%     alphas(alphas<0)=alphas(alphas<0)+360;
% end;
% while any(alphas>=360)
%     alphas(alphas>=360)=alphas(alphas>=360)-360;
% end;

iBeta=max(min(round(betas/betaStep)+1,nBeta),1);
beta=betaStep*(iBeta-1);
if ~isreal(beta)
    beta
    iBeta
    betaStep
    error('Complex beta');
end;
na=min(round(sind(beta)*nAlpha+.5),nAlpha);
q=alphas.*na/360-0.49;
iAlpha=max(min(round(q)+1,na),1);
% [beta na iAlpha]

ptrs=listIndices(sub2ind(size(listIndices),iAlpha,iBeta));
