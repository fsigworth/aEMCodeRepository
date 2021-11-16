function msk=hrpMaskDots(m);
% Make a 2D mask to delete the black dots in 20181216/sq05 data
% msk=1 where we want to blank the micrograph.
mf1=GaussFilt(GaussHP(m,.003),.05);

msk1=GaussFilt(mf1<-.11, .05)>.02;
% imags(msk1);
msk=(GaussFilt(msk1,.02)>.01);
% imags(msk2);

% imags(GaussFilt(m1u.*(1-msk2),.1));
