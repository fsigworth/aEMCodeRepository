function [means,medians,vals,rVals]=Radial3(w,org)
% function means = Radial3(w,org)
% function [means,medians,vals,rVals]=Radial3(w,org);
% Compute the circularly-averaged, radial component of the cubic volume w,
% with the origin taken to be org.  The number of points returned is
% nr=floor(min(size(w))/2).
% The assignment of radius values used to be ceil(r) but is now round(r). Thus
% means(1) = average of points in shell with radius in [0,.5), i.e. 0 only.
% means(2) = average of points in shell with radius in [.5, 1.5), etc.
% New addition fs July 2021:
% If desired, medians of shells are returned. Also, vals is a cell array, each
% element contains a vector of sorted values from each shell. (The median
% is taken as the middle value of that vector (or mean of the two middle
% values, if even.) rVals are the corresponding
% radius values.


sz = size(w);
if nargin<2 || numel(org)<1
    org=ceil((sz+1)/2);
end;
szmin=min(sz);
nr=floor(szmin/2);
R=Radius3(sz,org);
    [means, norm]=WeightedHisto(round(R)+1,w,nr);
    means=means./(norm+eps);

if nargout>1 % We're computing ordered entries
    rs=R(:);
    [rVals,rInds]=sort(rs);
    steps=find(diff(round(rVals)));
    steps=[0;steps];
    vals=cell(nr,1);
    medians=zeros(nr,1);
    for i=1:nr
        binVals=sort(w(rInds(steps(i)+1:steps(i+1))));
        nb=numel(binVals);
        if mod(nb,2)  % odd
            medians(i)=binVals((nb+1)/2);
        else
            lowCtr=nb/2;
            medians(i)=mean(binVals(lowCtr:lowCtr+1));
        end;
        vals{i}=binVals;
    end;
end;    
