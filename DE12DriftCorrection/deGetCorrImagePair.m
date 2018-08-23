function [im1 im2 navg]=deGetCorrImagePair(s,index,useMean)
% function [im1 imm]=deGetImagePair(s,index) 
% Given the DR camera movie s, return two copies of the image(index)
% corrected by the two sets of dark and bright references in s. If
% useMean=1 then im2 is instead returned as the mean of all the other images.
% 
% If index=0 or not given, then the returned im1 and im2
% are each the sum of all the images, corrected by different
% references.

if nargin<2
    index=0;
end;
if nargin<3
    useMean=0;
end;
nim=size(s.imgs,3);

cbr1=s.br1-s.dr1;
cbr2=s.br2-s.dr2;
if index>0
    im=single(s.imgs(:,:,index));
    im1=(im-s.dr1)./cbr1*mean(cbr1(:));
    if useMean
        imm=sum(single(s.imgs),3)-im;
        im2=(imm/(nim-1)-s.dr2)./cbr2*mean(cbr2(:));
        navg=nim-1;
    else
        im2=(im-s.dr2)./cbr2*mean(cbr2(:));
        navg=1;
    end;
else % index=0, so return the mean of all.
    im=sum(single(s.imgs),3);
    im1=(im-nim*s.dr1)./cbr1*mean(cbr1(:));
    im2=(im-nim*s.dr2)./cbr2*mean(cbr2(:));
end;

