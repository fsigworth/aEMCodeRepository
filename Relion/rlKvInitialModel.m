% rlKvInitialModel

outMap='KvRef128.mrc';
pixA=2.26;
boxSize=128;

m=arGetRefVolumes(pixA,boxSize,'KvMap.mat');
ShowSections(m);

WriteMRC(m,pixA,outMap);

