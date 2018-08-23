function [angleList, listIndices]=rsListHemisphereAngles(nAlpha, nBeta)
% Create an array t x 2 of angle values which are ~uniform on the upper
% hemisphere. The angles, in degrees, are returned as angleList(i,:)=[alpha
% beta], with alpha changing most rapidly. A lookup table of indices of
% alpha and beta steps, approximately an upper triangular matrix, is
% returned as listIndices(iAlpha,iBeta)=i.
% If nBeta is odd, we include a point at 90 degees; if even, not.
% e.g. nBeta = 2 yields 0 and 60 degrees.
% 
% Example:
% [angs,inds]=rsListHemisphereAngles(5,5);
% angs =
%             0            0
%             0         22.5
%           180         22.5
%             0           45
%            90           45
%           180           45
%           270           45
%             0         67.5
%            72         67.5
%           144         67.5
%           216         67.5
%           288         67.5
%             0           90
%            72           90
%           144           90
%           216           90
%           288           90
% inds =
%            1           2           4           8          13
%            0           3           5           9          14
%            0           0           6          10          15
%            0           0           7          11          16
%            0           0           0          12          17

corr=(mod(nBeta,2)+1)/2; % .5 even, 1 odd
betaStep=90/(nBeta-corr); % odd nBeta gives a final value at 90 degrees;
% even nBeta gives a final value of 90*(nBeta-1)/(nBeta-.5)), so that
% thetas will be evenly spaced when the lower hemisphere is included.
k=0;
angleList=zeros(nAlpha*nBeta,2);  % allocate an array bigger than we need.
listIndices=int32(zeros(nAlpha,nBeta));
for i=1:nBeta
  beta=(i-1)*betaStep;
  na=min(round(sind(beta)*nAlpha+.5),nAlpha);
  alphaStep=360/na;
%   [beta na]
  for j=1:na
      alpha=(j-1)*alphaStep;
      k=k+1;
      listIndices(j,i)=k;  % lookup
      angleList(k,:)=[alpha beta];
  end;
end;
angleList=angleList(1:k,:);  % truncate the list.

% % test code
% xs=sind(angleList(:,1)).*sind(angleList(:,2));
% ys=cosd(angleList(:,1)).*sind(angleList(:,2));
% plot(xs,ys,'.-');