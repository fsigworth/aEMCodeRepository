% rlxMakeScaledMicrographs
%
% [temporary fix] Copy the micrographs into the Merged/ directory, just scaling them
% to match the _v.mrc subtracted micrographs, but not padding.

% load Picking_9/allMis.mat

nmi=numel(allMis);

for i=1:nmi
    mi=allMis{i};
    m=ReadMRC([mi.imagePath mi.imageFilenames{1}]);
    mscl=(m-mi.imageMedian)*mi.imageNormScale;
    outName=[mi.procPath mi.baseFilename '_u.mrc'];
    disp([num2str(i) outName]);
    WriteMRC(mscl,mi.pixA,outName);
end
