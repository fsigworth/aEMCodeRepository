function [me1 me2 navgd]=deGetAlternateMeans(s)
% function [me1 me2 navgd]=deGetAlternateMeans(s)
% Similar to deGetImagePair, but returns the means of alternating frames of
% the movie in s.  The averages are corrected with references 1 and 2,
% respectively.  navgd is a 1x2 vector containing the numbers of frames
% averaged in each case.

nim=size(s.imgs,3);
ni1=ceil(nim/2);
ni2=nim-ni1;

cbr1=s.br1-s.dr1;  % corrected bright reference
cbr2=s.br2-s.dr2;

me1=(sum(single(s.imgs(:,:,1:2:nim)),3)-ni1*s.dr1)./cbr1*mean(cbr1(:))/ni1;
me2=(sum(single(s.imgs(:,:,2:2:nim)),3)-ni2*s.dr2)./cbr2*mean(cbr2(:))/ni2;    

navgd=[ni1 ni2];
