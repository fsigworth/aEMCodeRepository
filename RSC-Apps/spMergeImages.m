function spMergeImages(names)
if nargin<1
    names=[];
end;
pars.overridePixA=1.2156;
pars.overwrite=1;
pars.initialDefoci=1.5;
pars.removeOutliers=1;
MergeImages(names,pars);
